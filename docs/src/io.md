# `io` モジュール

TOML 設定ファイルの読み込みと CSV 出力を担当する。

---

## モジュール構成

```
io/
├── mod.rs     — サブモジュールの公開宣言
├── input.rs   — TOML 設定ファイルのデシリアライズ・検証
└── output.rs  — CSV 出力・サマリー表示
```

---

## `input.rs` — 設定ファイル

### `FlameConfig::from_file(path)`

TOML ファイルを読み込み、値の検証を行って `FlameConfig` を返す。
検証エラーは `anyhow::Error` として返る。

### 設定構造体

#### `[mechanism]`

| フィールド | 型 | デフォルト | 説明 |
|---|---|---|---|
| `file` | `String` | — | 機構ファイルのパス |
| `format` | `String` | `"cantera_yaml"` | ファイル形式（現在 `cantera_yaml` のみ有効） |

#### `[flame]`

未燃組成は **当量比モード** と **直接指定モード** のいずれか一方を使用する。
両方同時に指定するとバリデーションエラーになる。

| フィールド | 型 | デフォルト | 説明 |
|---|---|---|---|
| `pressure` | `f64` | `101325.0` | 圧力 [Pa] |
| `t_unburned` | `f64` | — | 未燃混合気温度 [K]（必須） |
| `domain_length` | `f64` | `0.02` | 計算領域長さ [m] |
| `fuel` | `{String → f64}` | — | 燃料モル分率マップ（当量比モード） |
| `oxidizer` | `{String → f64}` | — | 酸化剤モル分率マップ（当量比モード） |
| `equivalence_ratio` | `f64` | — | 当量比 φ（当量比モード） |
| `composition` | `{String → f64}` | — | 未燃混合気モル分率（直接指定モード） |

#### `[grid]`

| フィールド | 型 | デフォルト | 説明 |
|---|---|---|---|
| `initial_points` | `usize` | `20` | 初期格子点数（≥ 2） |
| `max_points` | `usize` | `500` | 最大格子点数（≥ initial_points） |
| `grad` | `f64` | `0.05` | GRAD 精細化基準（0 < grad ≤ 1） |
| `curv` | `f64` | `0.10` | CURV 精細化基準（0 < curv ≤ 1） |
| `ratio` | `f64` | `2.0` | 隣接セル幅の最大比 |

#### `[solver]`

| フィールド | 型 | デフォルト | 説明 |
|---|---|---|---|
| `atol` | `f64` | `1e-9` | Newton 絶対収束許容値 |
| `rtol` | `f64` | `1e-6` | Newton 相対収束許容値 |
| `max_newton_iter` | `usize` | `50` | Newton 最大反復回数 |
| `time_steps` | `usize` | `100` | 疑似過渡継続のステップ数 |
| `dt_initial` | `f64` | `1e-7` | 疑似過渡継続の初期タイムステップ [s] |
| `su_initial_guess` | `f64?` | — | 火炎速度初期推定 [m/s]（省略時 0.4 m/s） |
| `initial_profile` | `String?` | — | CSV 初期プロファイルのパス（省略時はシグモイド初期推定） |

#### `[output]`

| フィールド | 型 | デフォルト | 説明 |
|---|---|---|---|
| `file` | `String` | `"flame_solution.csv"` | 出力 CSV ファイルのパス |

### バリデーション規則

| 条件 | エラーメッセージに含む語 |
|---|---|
| `pressure > 0` | `pressure` |
| `t_unburned > 0` | `t_unburned` |
| `domain_length > 0` | `domain_length` |
| `composition` と `fuel/oxidizer/equivalence_ratio` の同時指定は不可 | — |
| 当量比モード選択時は `fuel`, `oxidizer`, `equivalence_ratio` 全て必須 | — |
| `equivalence_ratio > 0` | `equivalence_ratio` |
| `initial_points >= 2` | `initial_points` |
| `max_points >= initial_points` | `max_points` |
| `0 < grad ≤ 1` | `grad` |
| `0 < curv ≤ 1` | `curv` |
| `atol > 0` | `atol` |
| `rtol > 0` | `rtol` |

### TOML 例 — 当量比モード

```toml
[mechanism]
file = "data/h2o2.yaml"

[flame]
pressure = 101325.0
t_unburned = 300.0
domain_length = 0.02
equivalence_ratio = 1.0

[flame.fuel]
H2 = 1.0

[flame.oxidizer]
O2 = 0.21
N2 = 0.79

[grid]
initial_points = 20
max_points = 500
grad = 0.05
curv = 0.10
ratio = 2.0

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
time_steps = 100
dt_initial = 1.0e-7
su_initial_guess = 2.0        # 省略可: デフォルト 0.4 m/s
# initial_profile = "data/cantera_h2air_initial.csv"  # 省略時はシグモイド初期推定

[output]
file = "flame_solution.csv"
```

### TOML 例 — 直接組成指定モード

```toml
[flame]
t_unburned = 300.0
composition = { H2 = 0.3, O2 = 0.15, N2 = 0.55 }
```

### テスト (5 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_valid_composition_mode` | 直接指定モードが正常にパースされること |
| `test_valid_equivalence_ratio_mode` | 当量比モードが正常にパースされること |
| `test_negative_pressure_rejected` | `pressure < 0` でエラーになること |
| `test_zero_equivalence_ratio_rejected` | `equivalence_ratio = 0` でエラーになること |
| `test_composition_and_phi_together_rejected` | 両モード同時指定でエラーになること |
| `test_negative_domain_length_rejected` | `domain_length < 0` でエラーになること |

---

## `output.rs` — CSV 出力

### `write_csv(path, x, mech, grid, pressure)`

収束解ベクトル `x` を CSV ファイルに書き出す。

**列構成** (計 6 + 2K 列、K = 化学種数):

| 列名 | 単位 | 内容 |
|---|---|---|
| `z [m]` | m | 格子点座標 |
| `T [K]` | K | 温度 |
| `M [kg/m2/s]` | kg/(m²·s) | 質量流束（全行同一値） |
| `u [m/s]` | m/s | 流速 u = M / ρ |
| `rho [kg/m3]` | kg/m³ | 密度（理想気体則） |
| `hrr [W/m3]` | W/m³ | 熱発生率 HRR = −Σk ωk·hk |
| `X_<name>` × K | — | モル分率 Xk = (Yk/Wk) / Σ(Yj/Wj) |
| `Y_<name>` × K | — | 質量分率 Yk |

`M [kg/m2/s]` 列は `initial_guess_from_csv` が読み込む列と同一のフォーマットであるため、
本ソルバーの出力 CSV をそのまま次の計算（当量比スイープ等）の `initial_profile` として使用できる。

**数値フォーマット**: `z`, `M`, `u`, `rho`, `hrr`, `Xk`, `Yk` は `{:.6e}`、`T` は `{:.4}`。

### `print_summary(x, mech, grid, pressure)`

標準出力にサマリーを表示する。

```
---------------------------------------------------
  Laminar flame speed Su = X.XXXX m/s
  Max temperature       = XXXX.X K
  Grid points           = XXX
  Mass flux M           = X.XXXXe-X kg/(m²·s)
---------------------------------------------------
```

- `Su = M / ρ(j=0)` — 左境界の密度を使って質量流束から火炎速度に換算

### テスト (2 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_mole_fractions_pure_n2` | 純 N2 プロファイルで X_N2 = 1.0、ΣXk = 1.0、HRR ≈ 0 |
| `test_write_csv_column_count` | write_csv が 6 + 2K 列の CSV を生成すること |
