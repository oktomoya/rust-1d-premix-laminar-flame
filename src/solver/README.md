# `solver` モジュール

線形系の求解・Newton 法・疑似過渡継続を担当する。

---

## モジュール構成

```
solver/
├── mod.rs               — サブモジュールの公開宣言
├── banded.rs            — 帯行列 LU 分解・求解
├── dense.rs             — 密行列 LU 分解・求解 (faer)
├── newton.rs            — 減衰 Newton 法（直線探索付き）
└── pseudo_transient.rs  — 疑似過渡継続 (PT)
```

---

## `banded.rs` — 帯行列

### `BandedMatrix`

LAPACK 互換の一般帯行列 (GB) 形式で係数行列を保持する。

**格納レイアウト**: `data[col × ldab + (kl + ku + row − col)]`

- `ldab = 2×kl + ku + 1`（部分ピボッティングの fill-in 分の kl 行を余分に確保）

| メソッド | 説明 |
|---|---|
| `BandedMatrix::new(n, kl, ku)` | n×n 帯行列を確保（全ゼロ初期化） |
| `set(row, col, val)` | 要素を書き込む |
| `get(row, col)` | 要素を読む |
| `factor_in_place()` | 部分ピボッティング付き帯 LU 分解を in-place で実施。ピボット行インデックスを `ipiv` に格納。計算量 O(n·kl·ku) |
| `solve_factored(b)` | `factor_in_place` 済みの行列で A·x = b を解く（b を上書き）。前進代入 → 後退代入 |
| `solve(b)` | `factor_in_place` + `solve_factored` を一括実行 |

**fill-in の扱い**: ピボット行入れ替え後、行 k の非ゼロ範囲は `k+ku+kl` 列まで
拡張しうる。後退代入でも同じ幅を使用することで fill-in を正確に処理する。

### テスト (3 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_banded_tridiagonal_solve` | 5×5 三重対角系 (kl=ku=1)。厳密解 [5,8,9,8,5] |
| `test_banded_wider_band` | 4×4 非対称帯行列 (kl=2, ku=1)。厳密解 [1,1,1,1] |
| `test_banded_pivoting_required` | A[0][0]=0 でピボッティング必須な 3×3 系。厳密解 [1,2,2] |

---

## `dense.rs` — 密行列

### `DenseMatrix`

faer の部分ピボッティング LU で n×n 密行列を求解する。

帯行列では帯外に落ちる質量流束 M 列（解ベクトル末尾、全残差に結合）を
正確に扱うために使用する。

| メソッド | 説明 |
|---|---|
| `DenseMatrix::new(n)` | n×n ゼロ行列を確保 |
| `set(row, col, val)` | 要素を書き込む |
| `get(row, col)` | 要素を読む |
| `solve(b)` | faer の `partial_piv_lu().solve()` で A·x = b を解く（b を上書き） |

### テスト (2 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_dense_tridiagonal_solve` | 5×5 三重対角系。厳密解 [5,8,9,8,5] |
| `test_dense_with_global_coupling` | 最終列が全行に結合する 4×4 系（M 変数の模擬）。厳密解 [0.75,0.75,0.75,1.5] |

---

## `newton.rs` — 減衰 Newton 法

### アルゴリズム概要

```
for iter in 0..max_iter:
    評価: F(x)
    収束判定: ‖F‖ < atol·√n + rtol·‖x‖
    Jacobian: J = numerical_jacobian_dense(x, F)  ← 最大 max_jac_age 反復再利用
    線形系:   J·step = −F
    ステップクリップ: 物理スケール上限を超える成分を縮小
    直線探索: α を 1 から 0.5 刻みで最大 12 回試行
    更新:     x ← x + α·step
    物理クランプ: T, Y, M を有効範囲に制限
    失速検出: best ‖F‖ が 1% 以上改善しない反復が max_stall_iter 回続いたら受理
```

> **密 Jacobian を使う理由**: 帯行列近似では M 列（解ベクトル末尾）が帯外に落ち、
> 固有値制約の感度情報が失われる。`numerical_jacobian_dense` で全結合を捉えることで
> Newton ステップの精度が向上する。

### `NewtonConfig`

| フィールド | デフォルト | 説明 |
|---|---|---|
| `atol` | `1e-9` | 絶対収束許容値 |
| `rtol` | `1e-6` | 相対収束許容値 |
| `max_iter` | `50` | 最大反復回数 |
| `max_jac_age` | `3` | Jacobian の再利用上限反復数 |
| `max_stall_iter` | `5` | 失速と判定するまでの反復数 |

Jacobian の強制再計算: ‖F‖ が前回評価時から 10 倍以上改善した場合、
線形化点が大きく変わったとみなして `jac_age` をリセットする。

### ステップクリップ (`step_clip_factor`)

| 変数 | 最大変化量 |
|---|---|
| T (温度) | 500 K |
| Yk (質量分率) | 0.5 |
| M (質量流束) | 2 × max(|M|, 1) kg/(m²·s) |

### 物理クランプ (`clamp_physical`)

| 変数 | 有効範囲 |
|---|---|
| T | [1, 10000] K |
| Yk | [0, 1] |
| M | (0, ∞) （下限 1e-6） |

### テスト (2 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_newton_converges_immediately_at_zero_residual` | ‖F‖ ≈ 0 の一様 N2 プロファイルで iter 0 に収束すること |
| `test_newton_scaled_convergence_criterion` | ‖F‖ < atol·√n + rtol·‖x‖ が一様解で成立すること |

---

## `pseudo_transient.rs` — 疑似過渡継続

### 目的

Newton 法の収束域外（点火前の急勾配・スティッフ化学反応領域）から Newton 域へ
誘導する。TWOPNT スタイルの後退オイラー法を採用。

### アルゴリズム概要

```
for step in 0..n_steps:
    評価: F_steady(x)
    ‖F_steady‖ < 1e7 → Newton 域に達したので終了
    dt < 1e-13 → 収束断念して Newton へ

    rdt = 1/dt
    x_old = x

    内側 Newton (max_inner_iter 回):
        F_pt(x) = F_steady(x) + D/dt·(x − x_old)
        J_pt    = J_steady(x) + D/dt の対角（ρ·cp で T、ρ で Yk、1 で M）
        M 列をゼロクリア（M の疑似時間項は対角のみ）
        step = −J_pt⁻¹·F_pt
        ステップクリップ + 直線探索 + 物理クランプ

    外側受理: ‖F_steady(x_new)‖ ≤ 1.1×‖F_steady(x_old)‖
        成功 → dt *= dt_grow (最大 dt_max まで)
        失敗 → x = x_old, dt *= 0.5; 30 回連続失敗で中断
```

### `PseudoTransientConfig`

| フィールド | デフォルト | 説明 |
|---|---|---|
| `n_steps` | `300` | 最大外側ステップ数 |
| `dt_initial` | `1e-7` | 初期タイムステップ [s] |
| `dt_max` | `1e-3` | タイムステップ上限 [s] |
| `dt_grow` | `1.5` | 成功後のタイムステップ成長係数 |
| `max_inner_iter` | `5` | 内側 Newton の最大反復数 |

### 設計上の注意

**直線探索の初期化**: 内側直線探索の `best_norm` を ∞ で初期化することで、
J_pt が不定値（`D/dt < |λ_chem|`）のとき全 α で F_pt が増加しても
α = 0 に停滞しない。最小の有限 F_pt を与える α を返す。

**外側受理の 10% マージン**: 化学反応が非常に剛な点火フェーズでは、
良い内側 Newton ステップでも F_steady が一時的に増加することがある。
1.1 倍の許容マージンにより探索的なステップを通過させ、Newton 域への
到達を可能にする。

**D/dt 対角**: `ρ·cp` (T 方程式)、`ρ` (Yk 方程式)、`1` (M) で各変数の
次元が整合した疑似時間スケーリングを実現する（`residual.rs` の PT 埋め込みと同一）。

### テスト (2 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_pseudo_transient_reduces_residual` | T に小さな摂動を加えた N2 プロファイルに 5 ステップ適用し、T と M が物理的範囲内に収まること |
| `test_pseudo_transient_preserves_steady_solution` | 厳密定常解 (N2 均一) に 10 ステップ適用しても x が変化しないこと (max Δ < 1e-10) |
