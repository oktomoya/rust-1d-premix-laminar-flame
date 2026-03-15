# `transport` モジュール

1D 火炎ソルバーの輸送係数層。衝突積分・化学種個別輸送係数・混合則を担当する。

---

## モジュール構成

```
transport/
├── mod.rs                  — サブモジュールの公開宣言
├── collision_integrals.rs  — Neufeld (1972) 衝突積分多項式
├── species_props.rs        — 化学種単体の μk, λk, Dij
└── mixture.rs              — 混合則: μ_mix, λ_mix, Dkm
```

---

## `collision_integrals.rs` — 衝突積分

Chapman-Enskog 輸送理論で使用する換算衝突積分を Neufeld (1972) の多項式フィットで計算する。
引数 `t_star = T / (ε/kB)` は換算温度（無次元）。

### `omega22(t_star)` — Ω*(2,2)

粘性・熱伝導率の計算に使用。

```
Ω*(2,2) = A/T*^B + C/exp(D·T*) + E/exp(F·T*)
```

定数: A=1.16145, B=0.14874, C=0.52487, D=0.7732, E=2.16178, F=2.43787

### `omega11(t_star)` — Ω*(1,1)

二成分拡散係数の計算に使用。

```
Ω*(1,1) = A/T*^B + C/exp(D·T*) + E/exp(F·T*) + G/exp(H·T*)
```

定数: A=1.06036, B=0.15610, C=0.19300, D=0.47635, E=1.03587, F=1.52996, G=1.76474, H=3.89411

### 性質
- Ω*(2,2) > Ω*(1,1) (T* ≥ 1 で成立)
- 両積分ともに T* の単調減少関数
- T* → ∞ で剛体球極限 (~0.5–0.6) に収束
- 有効範囲: 0.3 ≤ T* ≤ 100（Neufeld 論文での検証範囲）

### テスト (7 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_omega22_polynomial_values` | 8点で Python 参照値と rtol=1e-10 |
| `test_omega11_polynomial_values` | 8点で Python 参照値と rtol=1e-10 |
| `test_omega22_vs_neufeld_table` | Neufeld Table I と rtol=0.2% |
| `test_omega11_vs_neufeld_table` | Neufeld Table I と rtol=0.2% |
| `test_monotone_decreasing` | T* 方向に単調減少 |
| `test_high_t_star_limit` | T*=100 で (0.4, 1.0) に収まる |
| `test_omega22_gt_omega11_high_t` | T* ≥ 1 で Ω22 > Ω11 |

---

## `species_props.rs` — 化学種輸送係数

### `viscosity(species, t)` → [Pa·s]

Chapman-Enskog 理論による純粋化学種の粘性。

```
μk = 2.6693e-6 × √(Wk[g/mol] × T) / (σk[Å]² × Ω*(2,2)(T*))
```

換算温度 `T* = T / (ε/kB)` は最小値 0.1 にクランプ。

### `thermal_conductivity(species, mu_k, cp_k, t, pressure)` → [W/(m·K)]

Mason-Monchick 式による純粋化学種の熱伝導率。並進モードと内部モードを分離して扱う。

| 種類 | 式 |
|---|---|
| 単原子 (Atom) | `λk = μk × 2.5 × cv_trans` |
| 線形・非線形 | `λk = μk × (2.5 × cv_trans + f_int × cv_int)` |

- `cv_trans = 3/2 × R/Wk`（並進比熱）
- `cv_int = cp_k - 5/2 × R/Wk`（内部比熱；回転＋振動）
- `f_int = ρk × D_kk / μk`（自己拡散比）

`f_int` は `binary_diffusion(k, k)` から計算する自己拡散係数を使用。簡略 Eucken 式（`f_int = 1.0`）に対して H₂O で約 13% の補正。

> **注意**: 旧来の簡略 Eucken 式 `μk × (cp_k + 1.25 R/Wk)` は `f_int = 1.0` に相当する。
> Mason-Monchick は `f_int ≈ 1.2–1.4` を使うため、特に極性分子（H₂O）で精度が向上する。

### `binary_diffusion(sp_i, sp_j, t, pressure)` → [m²/s]

BSL 式による二成分拡散係数。

```
Dij = 2.6280e-3 × T^1.5 / (P[atm] × σij[Å]² × Ω*(1,1)(T*) × √Wij[g/mol])
```

- `σij = (σi + σj) / 2`（算術平均）
- `εij = √(εi × εj)`（幾何平均）
- `Wij = 2WiWj/(Wi+Wj)`（換算分子量）

### テスト (7 テスト、全 PASS)

| テスト名 | 検証内容 |
|---|---|
| `test_viscosity_n2_300k` | μ_N2 300K vs Cantera 3.1.0（rtol=0.5%） |
| `test_viscosity_n2_1000k` | μ_N2 1000K vs Cantera 3.1.0（rtol=0.5%） |
| `test_binary_diffusion_h2n2_300k` | D_H2N2 300K vs Cantera 3.1.0（rtol=0.5%） |
| `test_binary_diffusion_h2n2_1000k` | D_H2N2 1000K vs Cantera 3.1.0（rtol=0.5%） |
| `test_thermal_conductivity_n2_300k` | λ_N2 300K vs NIST（rtol=10%） |
| `test_thermal_conductivity_h2_300k` | λ_H2 300K vs NIST（rtol=5%） |
| `test_thermal_conductivity_ar_300k` | λ_Ar 300K vs NIST（rtol=5%） |

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

混合平均拡散係数（各化学種 k に対して）。Cantera の `getMixDiffCoeffs` と一致する形式。

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
| `test_wilke_h2n2_300k` | μ_H2N2 (X=0.5/0.5) 300K（大 MW 差で Wilke 指数バグの回帰テスト、rtol=2%） |

---

## 数値的注意事項

- `T*` は最小値 0.1 にクランプ（超低温での発散を防ぐ）
- `sum_jk = 0` の場合（希薄モル分率など）、`Dkm = 0` を返す
- Wilke 則でモル分率が `< 1e-20` の成分はスキップ
- 二成分拡散係数の計算で `Dkj < 1e-30` の場合はクランプして除算を防ぐ

---

## 参考文献

- Neufeld, Janzen, Aziz (1972). *J. Chem. Phys.* 57, 1100.
- Bird, Stewart, Lightfoot (2002). *Transport Phenomena*, 2nd ed. Wiley. §9.3.
- Wilke (1950). *J. Chem. Phys.* 18, 517.
- Mason, Monchick (1962). *J. Chem. Phys.* 36, 1622.
