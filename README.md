# rust-1d-premix-laminar-flame

> **実験的プロジェクト。** Rust で 1D 予混合層流火炎ソルバーを実装する個人実験です。
> 動作する完成品ではなく、実装・検証の過程を記録することを目的としています。

燃料・酸化剤混合気の未燃条件（温度・圧力・当量比）を与えると、火炎全域の温度・
化学種質量分率の空間プロファイルおよび層流燃焼速度（固有値）を計算するソルバー。

---

## アプローチ

PREMIX / Cantera の手法に準拠：

- **支配方程式**: 1D 非一様格子上の化学種保存式 + エネルギー方程式
- **化学反応**: Arrhenius、Troe フォールオフ、PLOG 圧力依存反応を含む詳細反応機構
- **輸送**: 混合平均（Chapman-Enskog 粘性、二成分拡散、Wilke 混合則）
- **熱力学**: NASA 7 係数多項式
- **離散化**: 有限差分（対流：風上差分、拡散：中心差分）
- **ソルバー**: 減衰ニュートン法 + 擬似過渡継続法 + 適応格子精細化（GRAD/CURV）
- **入力**: Cantera YAML 反応機構ファイル、TOML 設定ファイル
- **出力**: z, T, u, ρ, Yk の CSV プロファイル

参考実装: Kee et al., "PREMIX" (Sandia, 1985) および Cantera `FreeFlame`。

---

## ディレクトリ構成

```
src/
├── chemistry/
│   ├── species.rs             # Species 構造体（NASA 係数、輸送パラメータ）
│   ├── mechanism.rs           # Mechanism: 化学種 + 反応
│   ├── thermo.rs              # NASA7 多項式: cp, h, s
│   ├── kinetics.rs            # Arrhenius、Troe フォールオフ、生成速度 ωk
│   └── parser/
│       ├── cantera_yaml.rs    # Cantera YAML パーサー  [実装済み]
│       └── chemkin.rs         # CHEMKIN-II パーサー    [スタブのみ]
├── transport/
│   ├── collision_integrals.rs # Neufeld (1972) Ω*(1,1), Ω*(2,2) 多項式フィット
│   ├── species_props.rs       # μk, λk, Dij（Chapman-Enskog）
│   └── mixture.rs             # μ_mix（Wilke）、λ_mix、Dkm（混合平均）
├── flame/
│   ├── state.rs               # 解ベクトルのレイアウト・インデックス
│   ├── residual.rs            # F(x): 化学種・エネルギー BVP + 境界条件
│   ├── refine.rs              # GRAD/CURV 適応格子精細化
│   └── solver_driver.rs       # 火炎ソルバー上位制御
├── solver/
│   ├── newton.rs              # 直線探索付き減衰ニュートン法
│   ├── pseudo_transient.rs    # 擬似過渡継続法
│   └── banded.rs              # 帯行列 LU 分解
└── io/
    ├── input.rs               # TOML 設定ファイル読み込み
    └── output.rs              # CSV 出力
```

---

## ビルド

```bash
cargo build
```

---

## 実行

```bash
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
T_unburned = 300.0
domain_length = 0.02

[grid]
initial_points = 20
max_points = 500
grad = 0.05
curv = 0.10

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
```

---

## テスト

全テストを実行:

```bash
cargo test
```

テスト名をモジュール別に一覧表示:

```bash
cargo test --lib -- --list
```

現在 38 テスト、全 PASS:

| モジュール | テスト数 | 検証内容 |
|---|---|---|
| `chemistry::parser::cantera_yaml` | 9 | 化学種・反応数、MW、Ea 単位、duplicate、falloff、三体、PLOG |
| `chemistry::thermo` | 12 | H2/O2/H2O/N2 の cp・h・s を 300K/1000K/2000K で Cantera 3.1.0 と比較（rtol=1e-8） |
| `chemistry::kinetics` | 4 | Kc（反応 11）、Troe 補正（9 ケース）、生成速度 ωk を Cantera 3.1.0 と比較 |
| `transport::collision_integrals` | 7 | Neufeld 多項式の算術精度、Table I クロス検証、単調減少性 |
| `transport::species_props` | 4 | μ_N2、D_H2N2 を 300K/1000K で Cantera 3.1.0 と比較（rtol=0.5%） |
| `transport::mixture` | 2 | Wilke μ_air（O2:N2=21:79）を 300K/1000K で Cantera 3.1.0 と比較（rtol=0.5%） |

すべての数値参照値は同一の `data/h2o2.yaml` 機構を用いた **Cantera 3.1.0** との
クロス検証によって確認済み。詳細は `VERIFICATION_REPORT.md` を参照。

---

## 参考文献

- Kee, Grcar, Smooke, Miller (1985). *PREMIX: A Fortran Program for Modeling Steady Laminar One-Dimensional Premixed Flames.* Sandia Report SAND85-8240.
- Kee, Rupley, Miller (1989). *Chemkin-II.* Sandia Report SAND89-8009.
- Neufeld, Janzen, Aziz (1972). *Empirical Equations to Calculate 16 of the Transport Collision Integrals.* J. Chem. Phys. 57, 1100.
- Cantera: https://cantera.org
