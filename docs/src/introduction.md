# rust-1d-premix-laminar-flame

> **実験的プロジェクト。** Rust で 1D 予混合層流火炎ソルバーを実装する個人実験です。

燃料・酸化剤混合気の未燃条件（温度・圧力・当量比）を与えると、火炎全域の温度・
化学種質量分率の空間プロファイルおよび層流燃焼速度（固有値）を計算するソルバー。

---

## アプローチ

PREMIX / Cantera の手法に準拠：

- **支配方程式**: 1D 非一様格子上の化学種保存式 + エネルギー方程式
- **化学反応**: Arrhenius、Troe フォールオフ、PLOG 圧力依存反応を含む詳細反応機構
- **輸送**: 混合平均輸送（Chapman-Enskog 粘性・拡散、Mason-Monchick 熱伝導率、Wilke 混合則）
- **衝突積分**: Monchick-Mason (1961) 2D テーブル（T* × δ*、極性分子対応）
- **熱力学**: NASA 7 係数多項式
- **離散化**: 有限差分（対流：風上差分、拡散：中心差分）
- **ソルバー**: 減衰ニュートン法（帯行列 LU + 枠付きシュア補完）+ 擬似過渡継続法 + 適応格子精細化（GRAD/CURV）
- **入力**: Cantera YAML 反応機構ファイル、TOML 設定ファイル
- **出力**: z, T, M, u, ρ, HRR, Xk, Yk の CSV プロファイル

参考実装: Kee et al., "PREMIX" (Sandia, 1985) および Cantera `FreeFlame`。

---

## 検証結果

H₂/air φ=1、300 K、1 atm（`data/h2o2.yaml`、Cantera 3.1.0 参照）:

| 量 | 本実装 | Cantera 3.1.0 | 誤差 |
|---|---|---|---|
| Su [m/s] | 2.3327 | 2.3354 | 0.1% |
| T_max [K] | 2357 | — | — |

残差 0.1% は風上差分の O(Δz) 打ち切り誤差によるもの

### H₂/air 当量比スイープ（`data/h2o2.yaml`、300 K、1 atm）

| φ | Su [m/s] | T_max [K] | 格子点数 |
|---|---|---|---|
| 0.8 | 1.649 | 2156 | 226 |
| 0.9 | 2.022 | 2284 | 221 |
| 1.0 | 2.333 | 2357 | 221 |
| 1.1 | 2.584 | 2364 | 222 |
| 1.2 | 2.787 | 2337 | 224 |

---

## ビルド・実行

```bash
cargo build
cargo run -- examples/hydrogen_air.toml
```

設定ファイル例（TOML）:

```toml
[mechanism]
file = "data/h2o2.yaml"
format = "cantera_yaml"

[flame]
pressure = 101325.0
fuel = { H2 = 1.0 }
oxidizer = { O2 = 0.21, N2 = 0.79 }
equivalence_ratio = 1.0
t_unburned = 300.0
domain_length = 0.02

[grid]
initial_points = 100
max_points = 1000
grad = 0.05
curv = 0.10

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
```

---

## 可視化

`scripts/plot_flame.html` — ブラウザで動作するインタラクティブな火炎プロファイルビューア（Plotly.js 使用）。

- **シングルモード**: 4 パネル（T+HRR、主要種、ラジカル種対数、流速+密度）
- **比較モード**: 最大 3 ファイルを重ねて比較（6 パネル + Su/T_max サマリーテーブル）
- 全トレースで `lines+markers` 表示（格子点位置を可視化）

---

## テスト

```bash
cargo test           # 全テスト（ユニット + 統合）
cargo test --lib     # ユニットテストのみ
```

現在 83 テスト、全 PASS（統合テスト 1 件含む）:

| モジュール | テスト数 | 検証内容 |
|---|---|---|
| `chemistry::parser::cantera_yaml` | 9 | 化学種・反応数、MW、Ea 単位、duplicate、falloff、三体、PLOG |
| `chemistry::parser::chemkin` | 13 | ELEMENTS/SPECIES パース、NASA7、輸送パラメータ、反応・falloff・三体効率 |
| `chemistry::thermo` | 12 | H2/O2/H2O/N2 の cp・h・s を 300K/1000K/2000K で Cantera 3.1.0 と比較（rtol=1e-8） |
| `chemistry::kinetics` | 4 | Kc、Troe 補正、生成速度 ωk を Cantera 3.1.0 と比較 |
| `transport::collision_integrals` | 11 | Neufeld 多項式、Monchick-Mason 2D テーブル、δ* 計算 |
| `transport::species_props` | 7 | μ_N2、D_H2N2、λ_N2/H2/Ar を Cantera 3.1.0 と比較（rtol ≤ 1%） |
| `transport::mixture` | 3 | Wilke μ_air（O2:N2=21:79）、μ_H2N2 を Cantera 3.1.0 と比較 |
| その他（flame、solver）| 23 | 残差ゼロ検証、Jacobian 正確性・帯行列 LU 正確性、擬似過渡継続、格子精細化など |
| `h2air_validation`（統合） | 1 | H₂/air φ=1: Su を Cantera 参照値 2.3354 m/s と比較（rtol=0.5%） |

---

## 参考文献

- Kee, Grcar, Smooke, Miller (1985). *PREMIX.* Sandia Report SAND85-8240.
- Kee, Rupley, Miller (1989). *Chemkin-II.* Sandia Report SAND89-8009.
- Monchick, Mason (1961). *Transport Properties of Polar Gases.* J. Chem. Phys. 35, 1676.
- Mason, Monchick (1962). *Heat Conductivity of Polyatomic and Polar Gases.* J. Chem. Phys. 36, 1622.
- Neufeld, Janzen, Aziz (1972). *Empirical Equations to Calculate 16 of the Transport Collision Integrals.* J. Chem. Phys. 57, 1100.
- Wilke (1950). *A Viscosity Equation for Gas Mixtures.* J. Chem. Phys. 18, 517.
- Cantera: https://cantera.org
