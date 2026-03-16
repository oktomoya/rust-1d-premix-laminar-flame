# `transport` モジュール

1D 火炎ソルバーの輸送係数層。衝突積分・化学種個別輸送係数・混合則を担当する。

---

## モジュール構成

```
transport/
├── mod.rs                  — サブモジュールの公開宣言
├── collision_integrals.rs  — Neufeld (1972) 多項式 + Monchick-Mason (1961) 2D テーブル
├── species_props.rs        — 化学種単体の μk, λk, Dij
└── mixture.rs              — 混合則: μ_mix, λ_mix, Dkm
```

---

## `collision_integrals.rs` — 衝突積分

Chapman-Enskog 輸送理論で使用する換算衝突積分を 2 つの方法で提供する。

### Neufeld (1972) 多項式フィット — 非極性化学種向け

引数 `t_star = T / (ε/kB)` は換算温度（無次元）。非極性化学種（δ* = 0）のみ有効。

#### `omega22(t_star)` — Ω*(2,2)

粘性・熱伝導率の計算に使用。

```
Ω*(2,2) = A/T*^B + C/exp(D·T*) + E/exp(F·T*)
```

定数: A=1.16145, B=0.14874, C=0.52487, D=0.7732, E=2.16178, F=2.43787

#### `omega11(t_star)` — Ω*(1,1)

二成分拡散係数の計算に使用。

```
Ω*(1,1) = A/T*^B + C/exp(D·T*) + E/exp(F·T*) + G/exp(H·T*)
```

定数: A=1.06036, B=0.15610, C=0.19300, D=0.47635, E=1.03587, F=1.52996, G=1.76474, H=3.89411

---

### Monchick-Mason (1961) 2D テーブル — 極性化学種対応

極性分子を含む化学種対に対し、換算双極子モーメント δ* の補正を施した衝突積分を
テーブル補間で計算する。Cantera の `MMCollisionInt` と同一データ・同一補間法。

テーブルは T* = 0.1–100（37 点）× δ* = 0–2.5（8 点）の 2D グリッド。

#### `omega22_mm(t_star, delta_star)` — Ω*(2,2)(T*, δ*)

δ* = 0 の場合はテーブルの第 0 列から補間（Neufeld 多項式と < 0.3% の差）。

#### `omega11_mm(t_star, delta_star)` — Ω*(1,1)(T*, δ*)

A*(T*, δ*) = Ω*(2,2) / Ω*(1,1) テーブルを使用して `omega22_mm / astar` で計算。

#### `delta_star_reduced(dipole_i, dipole_j, well_depth_ij, diameter_ij)` — δ*

```
δ*_ij = 0.5 × μi × μj / (4πε₀ × εij × σij³)
```

引数の単位: 双極子モーメント [Debye]、複合井戸深さ [K]、複合径 [Å]。

変換定数 3622.85 は単位変換を吸収する（1 D², 1 K, 1 Å³ での δ* の値）。

h2o2.yaml では H₂O のみ双極子モーメントを持つ（μ = 1.844 D → δ*_H₂O ≈ 1.22）。

#### 補間法

Cantera の `MMCollisionInt::quadInterp` と同一:
- δ* 方向: 隣接グリッド間の線形補間
- T* 方向: 3 点二次補間（log T* 座標）

---

### テスト (11 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_omega22_polynomial_values` | 8 点で Python 参照値と rtol=1e-10 |
| `test_omega11_polynomial_values` | 8 点で Python 参照値と rtol=1e-10 |
| `test_omega22_vs_neufeld_table` | Neufeld Table I と rtol=0.2% |
| `test_omega11_vs_neufeld_table` | Neufeld Table I と rtol=0.2% |
| `test_monotone_decreasing` | T* 方向に単調減少 |
| `test_high_t_star_limit` | T*=100 で (0.4, 1.0) に収まる |
| `test_omega22_gt_omega11_high_t` | T* ≥ 1 で Ω22 > Ω11 |
| `test_mm_omega22_delta0_matches_table` | δ*=0 でテーブル第 0 列と一致（rtol=1e-6） |
| `test_mm_omega11_delta0_matches_table` | δ*=0 で Neufeld 多項式と < 0.5% |
| `test_mm_omega22_h2o_1000k` | H₂O 1000K（δ*≈1.22）の範囲チェック |
| `test_delta_star_zero_for_nonpolar` | 非極性種で δ*=0 |

---

## `species_props.rs` — 化学種輸送係数

### `viscosity(species, t)` → [Pa·s]

Chapman-Enskog 理論による純粋化学種の粘性。

```
μk = 2.6693e-6 × √(Wk[g/mol] × T) / (σk[Å]² × Ω*(2,2)(T*, δ*k))
```

- δ*k > 0 の場合は `omega22_mm` を使用（H₂O のみ）
- δ*k = 0 の場合は Neufeld 多項式 `omega22` を使用

### `thermal_conductivity(species, mu_k, cp_k, t, pressure)` → [W/(m·K)]

完全な Mason-Monchick (1962) 式による熱伝導率。回転緩和補正（Zrot）を含む。

```
cv_rot:  0.0 (原子) | 1.0 (線形) | 1.5 (非線形)  [単位: R]
cv_int = cp/R - 2.5 - cv_rot   [振動モード、単位: R]
```

回転緩和係数 fz(T*):
```
fz(T*) = 1 + π^1.5/√T* · (0.5 + 1/T*)  + (π²/4 + 2)/T*
```

補正係数:
```
A  = 2.5 - f_int
B  = Zrot × fz(298/ε) / fz(T/ε)  +  (2/π) × (5/3 × cv_rot + f_int)
c1 = (2/π) × A / B
f_rot   = f_int × (1 + c1)
f_trans = 2.5 × (1 - c1 × cv_rot / 1.5)
λk = (μk/Wk) × R × (f_trans×1.5 + f_rot×cv_rot + f_int×cv_int)
```

### `binary_diffusion(sp_i, sp_j, t, pressure)` → [m²/s]

BSL 式による二成分拡散係数。

```
Dij = 2.6280e-3 × T^1.5 / (P[atm] × σij[Å]² × Ω*(1,1)(T*, δ*ij) × √Wij[g/mol])
```

- `σij = (σi + σj) / 2`（算術平均）
- `εij = √(εi × εj)`（幾何平均）
- `Wij = 2WiWj/(Wi+Wj)`（換算分子量）
- `δ*ij = delta_star_reduced(μi, μj, εij, σij)`

μi = 0 または μj = 0 の場合（H₂O と非極性種の対）は δ*ij = 0 となり、Neufeld 多項式を使用。

### テスト (7 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_viscosity_n2_300k` | μ_N2 300K vs Cantera 3.1.0（rtol=0.5%） |
| `test_viscosity_n2_1000k` | μ_N2 1000K vs Cantera 3.1.0（rtol=0.5%） |
| `test_binary_diffusion_h2n2_300k` | D_H2N2 300K vs Cantera 3.1.0（rtol=0.5%） |
| `test_binary_diffusion_h2n2_1000k` | D_H2N2 1000K vs Cantera 3.1.0（rtol=0.5%） |
| `test_thermal_conductivity_n2_300k` | λ_N2 300K vs Cantera 3.1.0（rtol=1%） |
| `test_thermal_conductivity_h2_300k` | λ_H2 300K vs Cantera 3.1.0（rtol=1%） |
| `test_thermal_conductivity_ar_300k` | λ_Ar 300K vs Cantera 3.1.0（rtol=1%） |

---

## `mixture.rs` — 混合則

### `mixture_viscosity(mech, mole_fractions, t)` → [Pa·s]

Wilke 混合則。

```
μ_mix = Σk Xk × μk / φk
```

相互作用係数:

```
φkj = [1 + (μk/μj)^0.5 × (Wj/Wk)^0.25]² / [√8 × √(1 + Wk/Wj)]
φk  = Σj Xj × φkj
```

### `mixture_thermal_conductivity(mech, mole_fractions, mass_fractions, t, pressure)` → [W/(m·K)]

算術平均・調和平均の相乗平均（Wassilijewa 混合則）。

```
λ_mix = 0.5 × (Σk Xk × λk  +  1 / Σk Xk/λk)
```

### `mixture_diffusion_coefficients(mech, mole_fractions, mass_fractions, t, pressure)` → `Vec<f64>` [m²/s]

混合平均拡散係数（各化学種 k に対して）。Cantera の `getMixDiffCoeffs` と一致。

```
Dkm = (1 - Yk) / Σ_{j≠k} (Xj / Dkj)
```

分子は Xk（モル分率）ではなく Yk（質量分率）を使用する。Cantera の定式化は
`(W̄ - Xk·Wk) / (W̄ · Σ Xj/Dkj) = (1 - Yk) / Σ Xj/Dkj`。

### `mass_to_mole_fractions(mech, y)` → `Vec<f64>`

質量分率 → モル分率の変換。`W̄ = 1 / Σk(Yk/Wk)` を使用。

### テスト (3 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_wilke_air_300k` | μ_air (O2:0.21, N2:0.79) 300K vs Cantera 3.1.0（rtol=0.5%） |
| `test_wilke_air_1000k` | μ_air (O2:0.21, N2:0.79) 1000K vs Cantera 3.1.0（rtol=0.5%） |
| `test_wilke_h2n2_300k` | μ_H2N2 (X=0.5/0.5) 300K（大 MW 差での Wilke 指数バグ回帰テスト、rtol=2%） |

---

## 数値的注意事項

- `T*` は最小値 0.1 にクランプ（超低温での発散を防ぐ）
- `sum_jk = 0` の場合（希薄モル分率など）、`Dkm = 0` を返す
- Wilke 則でモル分率が `< 1e-20` の成分はスキップ
- 二成分拡散係数の計算で `Dkj < 1e-30` の場合はクランプして除算を防ぐ

---

## 参考文献

- Neufeld, Janzen, Aziz (1972). *J. Chem. Phys.* 57, 1100.
- Monchick, Mason (1961). *J. Chem. Phys.* 35, 1676.
- Mason, Monchick (1962). *J. Chem. Phys.* 36, 1622.
- Bird, Stewart, Lightfoot (2002). *Transport Phenomena*, 2nd ed. Wiley. §9.3.
- Wilke (1950). *J. Chem. Phys.* 18, 517.
