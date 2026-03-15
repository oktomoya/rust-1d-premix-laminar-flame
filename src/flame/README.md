# `flame` モジュール

1D 火炎ソルバーのコア層。格子・解ベクトル・残差・ヤコビアン・格子精細化・
ソルバー制御を担当する。

---

## モジュール構成

```
flame/
├── mod.rs            — サブモジュールの公開宣言
├── domain.rs         — 非一様 1D 格子 (Grid)
├── state.rs          — 解ベクトルのレイアウトとインデックス
├── residual.rs       — 残差評価 F(x)
├── jacobian.rs       — 数値ヤコビアン
├── refine.rs         — GRAD/CURV 適応格子精細化
└── solver_driver.rs  — 火炎計算の上位制御
```

---

## `domain.rs` — 格子

### `Grid`

非一様 1D 格子。格子点座標を `Vec<f64>` で保持する。

| メソッド | 返り値 | 説明 |
|---|---|---|
| `Grid::uniform(length, n)` | `Grid` | 等間隔格子を生成 |
| `n_points()` | `usize` | 格子点数 J |
| `dz()` | `Vec<f64>` (長さ J-1) | セル幅 `dz[j] = z[j+1] - z[j]` |
| `z_mid()` | `Vec<f64>` (長さ J-1) | 中点座標 `(z[j] + z[j+1]) / 2` |
| `insert_midpoint(j)` | — | 区間 j に中点を挿入 |

---

## `state.rs` — 解ベクトルのレイアウト

### 解ベクトル構造

```
x = [T_0, Y_{0,0}, …, Y_{K-1,0},  T_1, Y_{0,1}, …,  …,  T_{J-1}, …, Y_{K-1,J-1},  M]
```

- 長さ: `(K+1) × J + 1`
- 各格子点に `NATJ = K+1` 変数（温度 1 + 化学種 K）
- 末尾の 1 要素は質量流束 M [kg/(m²·s)]（固有値）

| インデックス関数 | 返り値 |
|---|---|
| `natj(mech)` | K+1 |
| `solution_length(mech, nj)` | (K+1)×J + 1 |
| `idx_t(natj, j)` | 格子点 j の温度インデックス |
| `idx_y(natj, j, k)` | 格子点 j の化学種 k インデックス |
| `idx_m(natj, nj)` | 質量流束インデックス（末尾）|

### `FlameState<'a>` / `FlameStateMut<'a>`

解ベクトルへの高レベルアクセサ。

```rust
state.temperature(j)      // T[j]
state.species(k, j)       // Y_k[j]
state.mass_flux()         // M
state.y_slice(j)          // &Y[0..K] at point j
```

### `initial_guess`

シグモイドプロファイルによる初期推定値の生成。

```
T(z) = T_u + σ((z - z_c)/z_w) × (T_b - T_u)
Y(z) = Y_u + σ((z - z_c)/z_w) × (Y_b - Y_u)
```

### `initial_guess_from_csv(path, mech, grid)`

Cantera などの外部ソルバーが出力した CSV から初期推定値を読み込む。

- 必須列: `z_m`, `T_K`, `M_kg_m2s`, `Y_<name>` (各化学種)
- CSV の z 座標から現在の格子に線形補間し、Y は正規化（和を 1 に）する
- M は CSV 全点の平均値を使用
- `solver` 設定の `initial_profile` に CSV パスを指定すると疑似過渡継続（Phase 1）がスキップされ、Newton 法から直接開始できる

---

## `residual.rs` — 残差評価

### 支配方程式

**内部点 j = 1 … J-2:**

化学種方程式（k = 0 … K-2）:
```
F_yk = M × (Y_k,j - Y_k,j-1)/Δz_m
     + (jk_{j-1/2} - jk_{j-3/2}) / Δz_av
     - ωk × Wk
```

最終化学種は和の閉鎖式:
```
F_yK = ΣYk - 1
```

エネルギー方程式:
```
F_T = -M × cp × (T_j - T_j-1)/Δz_m
    + (λ_{j+1/2}×(T_{j+1}-T_j)/Δz_p - λ_{j-1/2}×(T_j-T_j-1)/Δz_m) / Δz_av
    - Σk (jk_{j-1/2} + jk_{j+1/2})/2 × cpk × (T_j - T_j-1)/Δz_m
    - Σk ωk × hk
```

**左境界 (j = 0):**
```
F_T  = T_0 - T_u
F_yk = M × (Y_k,0 - Y_k,u) + jk_{1/2}    (inlet flux BC)
```

**右境界 (j = J-1):**
```
F_T  = T_{J-1} - T_{J-2}    (zero gradient)
F_yk = Y_k,J-1 - Y_k,J-2
```

**固有値閉鎖:**
```
F_M = T(j_fix) - T_fix
```

### `FlameConfig`

| フィールド | 単位 | 内容 |
|---|---|---|
| `pressure` | Pa | 圧力 |
| `t_unburned` | K | 未燃温度（左境界） |
| `y_unburned` | — | 未燃質量分率 |
| `z_fix` | m | 温度固定点の位置 |
| `t_fix` | K | 固定温度（固有値閉鎖） |

### `eval_residual(x, rhs, mech, grid, config, x_old, rdt)`

- 定常計算: `x_old = None`, `rdt = 0.0`
- 疑似過渡継続: `x_old` と `rdt = 1/Δt` を渡すと内部点に backward-Euler 項を加算:
  ```
  T 方程式:  + rdt × ρ × cp × (T - T_old)    [W/m³]
  Yk 方程式: + rdt × ρ × (Yk - Yk_old)        [kg/(m³·s)]
  ```
  密度・熱容量のスケーリングにより各方程式が次元整合する。

### 拡散フラックス

混合平均拡散係数による拡散フラックスに補正速度を適用して ΣJk = 0 を強制:
```
jk_raw = -ρ × Dkm × dYk/dz
jk     = jk_raw - Yk × Σj jk_raw,j     (upwind Yk, 左セル値)
```

### 単位系

- 濃度は mol/cm³ で `production_rates` に渡す（CGS、A値と整合）
- ωk の出力 mol/(cm³·s) × 1e6 → mol/(m³·s)、さらに × Wk で kg/(m³·s)

### テスト (2 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_residual_zero_for_uniform_n2` | 均一 N2 プロファイルで \|F\| < 1e-8 |
| `test_pseudo_transient_zero_when_x_eq_x_old` | x_old = x のとき PT 項がゼロ |

---

## `jacobian.rs` — 数値ヤコビアン

前進差分による数値ヤコビアン `J = ∂F/∂x` を帯行列形式で返す。

- 帯域幅: `kl = ku = 2×NATJ`（隣接ブロック間結合 + M 列の影響を吸収）
- 総対角数: `2×NATJ + 1`
- 摂動幅: `h = √(2ε_machine) × max(|x_j|, 1)`
- LAPACK `dgbsv` 互換の列優先帯行列 (`BandedMatrix`) に格納

---

## `refine.rs` — 適応格子精細化

### `RefineCriteria`

| フィールド | デフォルト | 説明 |
|---|---|---|
| `grad` | 0.05 | 勾配基準（0–1; 小さいほど細かい） |
| `curv` | 0.10 | 曲率基準（0–1; 小さいほど細かい） |
| `ratio` | 2.0 | 隣接セル幅の最大比（超えたら挿入） |
| `max_points` | 500 | 最大格子点数 |

TOML で `[grid]` セクションに記述することで上書き可能。

### 精細化判定 (Cantera `refine.cpp` 準拠)

区間 j に中点を挿入する条件（T と全化学種 Yk について評価）:

```
GRAD:  |φ_{j+1} - φ_j|          > grad × (φ_max - φ_min)
CURV:  |slope[j+1] - slope[j]|  > curv × (slope_max - slope_min)
RATIO: dz[j] > ratio × dz[j-1]  または  dz[j-1] > ratio × dz[j]
```

- **GRAD**: 値の変化量がドメイン全体のレンジに対して大きすぎる区間を検出。
- **CURV**: 一階微分の変化量（`slope[j] = dφ/dz`）がドメイン全体の傾き範囲に対して大きすぎる区間を検出。両辺とも単位は φ/m で次元整合。
- **RATIO**: 隣接セル幅の比が `ratio` を超える場合に挿入。炎帯で局所的に細かくなった格子が隣の粗いセルへ連鎖的に伝播するため、滑らかな格子遷移を保証する。


### `refine(x, mech, grid, criteria)` → `Option<(Grid, Vec<f64>)>`

精細化が必要な点がなければ `None`。挿入は逆順で行いインデックスの安定性を保つ。
新しい格子点の解は線形補間で生成する。質量流束 M は末尾要素として挿入後に再付加する。

---

## `solver_driver.rs` — ソルバー制御

### `run_flame(config)` の処理フロー

```
1. 機構ファイルの読み込み (Cantera YAML または CHEMKIN-II)
2. 未燃・既燃組成の計算（当量比 φ、または直接指定）
3. 断熱火炎温度 T_ad の推定（エンタルピー収支）
4. 初期格子の生成（均一、initial_points 点）
5. Phase 1: 疑似過渡継続 (n_steps, dt_initial)
     ※ initial_profile が指定されている場合はスキップ
       → Cantera 収束解を初期値とする場合などに使用
6. Phase 2: Newton 法 + 適応格子精細化のループ（最大 20 パス）
   a. Newton 法で収束（scaled residual < rtol）
   b. GRAD/CURV/RATIO 基準で格子精細化
   c. 精細化不要 (refine が None) になるまで繰り返す
   d. 精細化ごとに z_fix を T プロファイルから再計算
7. CSV 出力（z, T, u, ρ, HRR, X_*, Y_* 列）
```

### 格子精細化の収束挙動

H2/air φ=1 の典型的な実行例（Cantera CSV 初期値、grad=0.05, curv=0.10, ratio=2.0）:

| pass | 格子点数 | ‖F‖ 収束値 |
|---|---|---|
| 0 | 100 | 6e-4 |
| 1 | 111 | 1e-3 |
| … | … | … |
| 8 | 221 | 9e-3 |

最終 Su = 2.3313 m/s（Cantera 参照値 2.3354 m/s、誤差 0.2%）。

### 未実装 (TODO)

| 関数 | 状況 |
|---|---|
| `estimate_burned_composition` | 未燃組成をそのまま返すプレースホルダー（完全燃焼の元素バランスが必要） |
| `compute_compositions` | 当量比を使った厳密な燃料・酸化剤混合計算が未対応 |
| `z_fix` の動的追跡 | **実装済み** — `find_z_fix` が各 Newton パス後に T プロファイルから再計算 |
