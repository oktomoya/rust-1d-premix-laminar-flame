# `chemistry` モジュール

1D 火炎ソルバーの化学反応層。化学種データの保持・熱力学計算・反応速度計算・
機構ファイルのパースを担当する。

---

## モジュール構成

```
chemistry/
├── mod.rs            — サブモジュールの公開宣言
├── species.rs        — 化学種データ構造 (Species, NasaPoly, TransportParams)
├── mechanism.rs      — 機構全体のデータ構造 (Mechanism, Reaction, RateType)
├── thermo.rs         — NASA7 多項式による熱力学関数
├── kinetics.rs       — 反応速度計算・生成速度
└── parser/
    ├── mod.rs
    ├── cantera_yaml.rs — Cantera YAML パーサー ✅ 実装済み
    └── chemkin.rs      — CHEMKIN-II パーサー   ✅ 実装済み
```

---

## `species.rs` — 化学種データ構造

### `Species`
1 化学種が持つ全データをまとめた構造体。

| フィールド | 型 | 単位 | 内容 |
|---|---|---|---|
| `name` | `String` | — | 化学種名 (例: `"H2"`) |
| `molecular_weight` | `f64` | kg/mol | 分子量 |
| `composition` | `HashMap<String, f64>` | — | 元素組成 (例: `{"H": 2.0, "O": 1.0}`)。キーは大文字正規化済み |
| `nasa_low` | `NasaPoly` | — | 低温域 NASA7 多項式 |
| `nasa_high` | `NasaPoly` | — | 高温域 NASA7 多項式 |
| `transport` | `TransportParams` | — | 輸送パラメータ |

`nasa_poly(t)` メソッドが温度 `t` に応じて `nasa_low` / `nasa_high` を自動選択する
(切り替え温度 = `nasa_low.t_high` = T_mid)。

### `NasaPoly`
NASA 7係数多項式 1レンジ分。有効温度範囲 `[t_low, t_high]` と係数配列 `coeffs[7]`。

```
cp/R = a[0] + a[1]·T + a[2]·T² + a[3]·T³ + a[4]·T⁴
h/RT = a[0] + a[1]/2·T + a[2]/3·T² + a[3]/4·T³ + a[4]/5·T⁴ + a[5]/T
s/R  = a[0]·ln(T) + a[1]·T + a[2]/2·T² + a[3]/3·T³ + a[4]/4·T⁴ + a[6]
```

### `TransportParams`
Lennard-Jones パラメータと分子形状 (輸送係数モジュールで使用)。

| フィールド | 単位 |
|---|---|
| `geometry` | `Atom` / `Linear` / `Nonlinear` |
| `well_depth` | K (ε/k_B) |
| `diameter` | Å (σ) |
| `dipole_moment` | Debye |
| `polarizability` | Å³ |
| `rot_relax` | — (298 K における回転緩和衝突数 Z_rot) |

---

## `mechanism.rs` — 機構データ構造

### `RateType`
反応速度定数の種別を enum で表現する。

```rust
pub enum RateType {
    Arrhenius { a: f64, b: f64, ea: f64 },
    Falloff { high: Arrhenius, low: Arrhenius, troe: Option<TroeParams> },
    Plog { rates: Vec<(f64, Arrhenius)> },
}
```

- **`Arrhenius`**: 通常の Arrhenius 式。三体反応も含む (`third_body` フィールドが `Some` のとき)。
  - `ea` の単位は **J/mol** (パーサーが変換済み)
  - `a` の単位はパース元ファイルのまま (Cantera YAML は CGS: cm³/mol/s 等)
- **`Falloff`**: 高圧限界 / 低圧限界の Arrhenius に加え、オプションの Troe 補正係数。
  Lindemann 形式は `troe = None`。
- **`Plog`**: 圧力依存反応。`(圧力 [Pa], Arrhenius)` のリストで定義される。

### `TroeParams`
Troe 補正のパラメータ。Cantera YAML の `Troe:` ブロックに対応。

```
Fcent = (1-A)·exp(-T/T3) + A·exp(-T/T1) + exp(-T2/T)
```

フィールド: `a`, `t3`, `t1`, `t2: Option<f64>` (T2 は省略可)。

### `ThirdBodySpec`
三体衝突パートナーの拡張効率を保持する。

```rust
pub struct ThirdBodySpec {
    pub efficiencies: Vec<(usize, f64)>,  // (species_index, eff)
}
```

デフォルト効率は 1.0 (暗黙)。非デフォルト種のみ格納する。
`[M]_eff = Σ[Ci] + Σ(eff_k - 1)·[Ck]` で有効第三体濃度を計算する。

### `Reaction`

| フィールド | 型 | 内容 |
|---|---|---|
| `reactants` | `Vec<(usize, f64)>` | (種インデックス, 化学量論係数) |
| `products` | `Vec<(usize, f64)>` | 同上 |
| `rate` | `RateType` | 速度定数の種別 |
| `reversible` | `bool` | `<=>` なら `true` |
| `third_body` | `Option<ThirdBodySpec>` | 三体効率 |
| `duplicate` | `bool` | `duplicate: true` フラグ |

---

## `thermo.rs` — NASA7 熱力学関数

`R_UNIVERSAL = 8.314462618 J/(mol·K)` を定数として定義。

| 関数 | 出力単位 | 説明 |
|---|---|---|
| `cp_over_r(species, t)` | 無次元 | cp/R |
| `cp_species(species, t)` | J/(kg·K) | 定圧比熱 |
| `enthalpy_species(species, t)` | J/kg | 比エンタルピー |
| `enthalpy_molar(species, t)` | J/mol | モルエンタルピー |
| `entropy_species(species, t)` | J/(kg·K) | 標準エントロピー (1 atm 基準) |
| `cp_mixture(species, y, t)` | J/(kg·K) | 質量分率 y[k] で混合 cp |
| `mean_molecular_weight(species, y)` | kg/mol | 混合平均分子量 (1/Σ(Yk/Wk)) |
| `density(p, t, mean_mw)` | kg/m³ | 理想気体則による密度 |

> **注意**: `entropy_species` は J/(kg·K) を返す。`kinetics.rs::equilibrium_constant`
> では `sp.molecular_weight` を乗じて J/(mol·K) に戻してから使用している。

---

## `kinetics.rs` — 反応速度・生成速度

### `arrhenius(a, b, ea, t)`
```
k(T) = A · T^b · exp(−Ea / (R·T))
```
`ea` は J/mol。

### `production_rates(mech, t, concentrations, pressure)`
全化学種のモル生成速度 ωk [mol/(cm³·s)] を返す。`concentrations[k] = ρ·Yk/Wk·1e-6` [mol/cm³]。
A 値が CGS 単位 (cm³/mol/s) で格納されているため、濃度も mol/cm³ で渡す必要がある。

各反応について以下の順序で処理する:

1. **速度定数 kf** を `RateType` に応じて計算:
   - `Arrhenius` → `arrhenius(a, b, ea, t)`
   - `Falloff` → `falloff_rate(...)` (後述)
   - `Plog` → `plog_rate(...)` (後述)

2. **三体補正**:
   - `Falloff` の場合: `falloff_rate` 内で既に [M] を取り込み済みのため **追加乗算なし**
   - それ以外で `third_body.is_some()` の場合: `kf × [M]_eff`

3. **平衡定数 Kc** → 逆反応速度定数 `kr = kf / Kc`

4. **進行速度** `q = kf·Π[Ri]^ν_i − kr·Π[Pi]^ν_i`

5. **生成速度** `ωk += Σ ν_prod · q − Σ ν_react · q`

### `falloff_rate` — Troe/Lindemann フォールオフ

```
k∞ = Arrhenius(high, T)
k0 = Arrhenius(low,  T)
Pr = k0·[M] / k∞
F  = troe_broadening(Troe, T, Pr)   ← Troe あり
   = 1.0                             ← Lindemann (Troe なし)
kf = k∞ · Pr/(1+Pr) · F
```

Troe ブロードニング係数 F の計算:
```
Fcent = (1−A)·exp(−T/T3) + A·exp(−T/T1) + exp(−T2/T)
c = −0.4 − 0.67·log10(Fcent)
n =  0.75 − 1.27·log10(Fcent)
f1 = (log10(Pr) + c) / (n − 0.14·(log10(Pr) + c))
F  = 10^(log10(Fcent) / (1 + f1²))
```
Cantera の `TroeRate::updateTemp` と同一の式。

### `plog_rate` — PLOG 圧力依存反応

圧力方向に対数線形補間:
```
log k = log k1 + (log P − log P1) · (log k2 − log k1) / (log P2 − log P1)
```
圧力範囲外はクランプ (最低 / 最高圧エントリを使用)。

### `equilibrium_constant(mech, rxn_idx, t)`
ΔG°/RT から Kp を求め、モル数変化 Δν を補正して Kc に変換する。
```
Kp = exp(−ΔG°/RT) = exp(−Σ νk·(hk − T·sk) / RT)
Kc = Kp · (P_atm / RT)^Δν
```

---

## `parser/cantera_yaml.rs` — Cantera YAML パーサー ✅

`parse_file(path)` または `parse_str(content)` で `Mechanism` を返す。

### パース手順
1. `units:` ブロックから活性化エネルギー単位を読む (デフォルト: `cal/mol`)
2. `species:` セクション → 各種 MW・NASA7・輸送パラメータ を構築
3. `reactions:` セクション → 反応タイプ判定 → `RateType` を構築

### 活性化エネルギー単位変換 (→ J/mol)

| YAML 値 | 変換係数 |
|---|---|
| `cal/mol` | × 4.184 |
| `kcal/mol` | × 4184.0 |
| `J/mol` | × 1.0 |
| `kJ/mol` | × 1000.0 |
| `K` | × 8.314462618 (= R) |
| `J/kmol` | × 1e-3 |

### 反応タイプ判定ロジック

```
type: pressure-dependent-Arrhenius  →  Plog
rate-constant: が null              →  Falloff  (high-P / low-P / Troe を読む)
それ以外                             →  Arrhenius
```

三体反応の判定: `type: three-body` または式中に `(+M)` / `+ M` を含む場合、
`efficiencies:` マップから `ThirdBodySpec` を構築する。

### 化学量論係数パース
- `2 OH` (スペース区切り係数)
- `2OH` (係数が名前に直接付く)
- 係数なし → 1.0
- `(+M)` / `+ M` / 孤立した `M` トークンはスキップ

### 分子量計算
`composition:` マップのキーを大文字正規化してテーブル参照。
対応元素: H, O, N, C, AR, HE, S, CL, F, BR。

### テスト (9 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_parse_h2o2_species_count` | 10 種 |
| `test_parse_h2o2_reaction_count` | 29 反応 |
| `test_molecular_weights` | H2, O2, AR, N2 の MW |
| `test_activation_energy_units_cal_mol` | 反応3 Ea = 6260 cal/mol → 26191.84 J/mol |
| `test_duplicate_flag` | 反応 24–29 に duplicate=true |
| `test_falloff_reaction_22` | high.A, low.A, Troe.A |
| `test_three_body_reaction_1` | 第三体あり・効率非空 |
| `test_case_insensitive_elements` | 小文字元素名 (`{h:2, o:1}`) |
| `test_plog_parsing` | 2 エントリの PLOG 反応 |

Cantera 3.1.0 Python API との数値クロスバリデーション (`scripts/cross_validate.py`):
164 チェック全 PASS (MW・NASA7係数・全29反応の A/b/Ea・Troe 4パラメータ・三体効率)。

---

## `parser/chemkin.rs` — CHEMKIN-II パーサー ✅

### 対応ファイル形式

| ファイル | 内容 |
|---|---|
| `mech.inp` | ELEMENTS / SPECIES / REACTIONS ブロック |
| `therm.dat` | 固定列 4 行形式 NASA7 多項式 (GRI-Mech 準拠) |
| `tran.dat` | 1 行 1 種の輸送パラメータ |

### エントリポイント

```rust
// 文字列から直接パース
parse_mechanism(mech_src: &str, therm_src: Option<&str>, tran_src: Option<&str>) -> Result<Mechanism>

// ファイルパスから読み込み
parse_file(mech_path: &str, therm_path: Option<&str>, tran_path: Option<&str>) -> Result<Mechanism>
```

- `therm_src` / `therm_path` が `None` の場合は `mech.inp` 内の `THERMO` ブロックを参照
- `tran_src` / `tran_path` が `None` の場合は輸送データなし（デフォルト値）

### パース手順

1. **ELEMENTS ブロック**: `ELEMENTS … END` を抽出、大文字正規化
2. **SPECIES ブロック**: `SPECIES … END` を抽出、宣言順を保持
3. **therm.dat**: `1`/`2`/`3`/`4` 行インジケータ（列 79）で 4 行セットを認識し、固定列から NASA7 係数・温度範囲・元素組成・MW を読む
4. **tran.dat**: `名前 geometry ε/k σ μ α Z_rot` の空白区切り形式
5. **REACTIONS ブロック**: 反応式行 + 修飾子行（`LOW/`, `TROE/`, `REV/`, 効率行, `DUPLICATE`）をレコードに集約して変換

### 反応修飾子

| 修飾子 | 形式 | 内容 |
|---|---|---|
| `LOW/ A b Ea /` | falloff 低圧限界 Arrhenius |
| `TROE/ a T3 T1 [T2] /` | Troe ブロードニング (T2 省略可) |
| `REV/ A b Ea /` | 明示的逆反応（構造に保持、現在は未使用） |
| `Species/ eff /` | 三体衝突効率 (複数まとめて 1 行可) |
| `DUPLICATE` または `DUP` | 重複反応フラグ |

- `(+M)` または `+ M` を含む方程式は自動的に三体反応と判定
- `LOW/` があれば `RateType::Falloff`、なければ `RateType::Arrhenius`
- Ea はデフォルト cal/mol として J/mol に変換 (× 4.184)
- PLOG 形式は CHEMKIN-II 非標準のため非対応（Cantera YAML パーサーを使用）

### テスト (13 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_parse_elements` | H, O, AR 等の元素抽出 |
| `test_parse_species_list` | 9 種、宣言順 |
| `test_thermo_h2_coefficients` | H2 低温・高温の a[0] と温度範囲 |
| `test_thermo_molecular_weights` | H2, O2, H2O, OH, H, O, AR の MW |
| `test_thermo_missing_species_error` | 未知種でエラー |
| `test_transport_geometry` | H=Atom, H2=Linear, H2O=Nonlinear |
| `test_transport_values` | N2 の ε/k, σ, Z_rot |
| `test_reaction_count` | 5 反応 |
| `test_arrhenius_reaction` | H + O2 <=> O + OH の A, b, Ea |
| `test_falloff_reaction` | H + O2 (+M) の high.A, low.A, TROE.a |
| `test_third_body_efficiencies` | H2O 効率 = 10.6 |
| `test_three_body_reaction` | 2H + M <=> H2 + M の A 値 |
| `test_parse_file_roundtrip` | `data/chemkin_test/` 3 ファイル読み込み、9 種・5 反応・H2O MW・形状 |

テストデータ: `data/chemkin_test/{mech.inp, therm.dat, tran.dat}`
