# Getting Started

## 必要なもの

- [Rust](https://rustup.rs/) 1.75 以上（`rustup` でインストール推奨）
- Git

確認コマンド:
```bash
rustc --version   # rustc 1.75.0 以上
cargo --version
```

---

## インストール

```bash
git clone https://github.com/oktomoya/rust-1d-premix-laminar-flame.git
cd rust-1d-premix-laminar-flame
cargo build --release
```

初回ビルドは依存クレートのコンパイルで数分かかります。

---

## 最初の計算を走らせる

同梱の H₂/air φ=1.0 設定ファイルで計算してみましょう:

```bash
cargo run --release -- examples/hydrogen_air.toml
```

正常に収束すると、こんなサマリーが表示されます:

```
---------------------------------------------------
  Laminar flame speed Su = 2.3327 m/s
  Max temperature       = 2356.9 K
  Grid points           = 221
  Mass flux M           = 2.8294e-1 kg/(m²·s)
---------------------------------------------------
Solution written to /tmp/h2air_phi_1.0.csv
```

---

## 出力ファイルを可視化する

`scripts/plot_flame.html` をブラウザで開き、生成された CSV ファイルをドラッグ＆ドロップします。

または Python スクリプト:

```bash
pip install pandas matplotlib numpy
python scripts/plot_flame.py /tmp/h2air_phi_1.0.csv
```

---

## 設定ファイルの書き方

`examples/hydrogen_air.toml` を複製して編集するのが最も簡単です。

```toml
[mechanism]
file = "data/h2o2.yaml"       # Cantera YAML 反応機構ファイル
format = "cantera_yaml"

[flame]
pressure = 101325.0           # 圧力 [Pa]
fuel     = { H2 = 1.0 }       # 燃料組成（モル分率）
oxidizer = { O2 = 0.21, N2 = 0.79 }
equivalence_ratio = 1.0       # 当量比 φ
t_unburned = 300.0            # 未燃温度 [K]
domain_length = 0.02          # 計算領域長さ [m]

[grid]
initial_points = 100          # 初期格子点数
max_points = 1000             # 最大格子点数
grad = 0.05                   # GRAD 精細化基準（小さいほど細かい）
curv = 0.10                   # CURV 精細化基準

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
su_initial_guess = 2.3        # 火炎速度初期推定 [m/s]（省略可）

[output]
file = "/tmp/h2air_phi_1.0.csv"
```

### 当量比を変えるには

`equivalence_ratio` を変えるだけです。収束が遅い場合は、隣の φ で計算した CSV を
`initial_profile` に指定すると擬似過渡継続フェーズをスキップして高速に収束します:

```toml
[solver]
su_initial_guess = 1.6
initial_profile = "/tmp/h2air_phi_0.9.csv"   # 前の計算結果を初期値に使用
```

---

## テストを実行する

```bash
cargo test                    # 全 83 テスト（ユニット + 統合）
cargo test --lib              # ユニットテストのみ（統合テストを除く）
cargo test -- --nocapture     # 標準出力を表示しながら実行
```

統合テスト (`tests/h2air_validation.rs`) はフルソルバーを走らせるため、
数秒かかります。

---

## よくある問題

### `cargo run` が遅い

`--release` フラグを忘れていませんか？デバッグビルドは最適化なしなので
化学反応計算が 10 倍以上遅くなります。

### Newton が収束しない

- `su_initial_guess` を実際の火炎速度に近い値にしてみてください
- `initial_profile` に収束済み CSV を指定して PT フェーズをスキップしてみてください
- `grad` / `curv` を少し大きく（0.1 / 0.2）して格子を粗くしてみてください

### 出力 CSV が空 / 途中で止まる

ターミナルの出力に `Newton: stalled` が多発している場合は PT フェーズで詰まっています。
`time_steps` を増やすか（例: 1000）、`dt_initial` を小さく（例: `1e-8`）してみてください。
