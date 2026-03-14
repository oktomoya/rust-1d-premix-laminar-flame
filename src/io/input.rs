use anyhow::Result;
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Debug, Deserialize)]
pub struct FlameConfig {
    pub mechanism: MechanismConfig,
    pub flame: FlameSectionConfig,
    pub grid: GridConfig,
    pub solver: SolverConfig,
    pub output: OutputConfig,
}

#[derive(Debug, Deserialize)]
pub struct MechanismConfig {
    pub file: String,
    #[serde(default = "default_format")]
    pub format: String,
}

fn default_format() -> String { "cantera_yaml".to_string() }

#[derive(Debug, Deserialize)]
pub struct FlameSectionConfig {
    #[serde(default = "default_pressure")]
    pub pressure: f64,
    /// Equivalence-ratio mode: specify fuel, oxidizer, and equivalence_ratio.
    pub fuel: Option<HashMap<String, f64>>,
    pub oxidizer: Option<HashMap<String, f64>>,
    pub equivalence_ratio: Option<f64>,
    /// Direct composition mode: specify mole fractions directly (bypasses φ).
    pub composition: Option<HashMap<String, f64>>,
    pub t_unburned: f64,
    #[serde(default = "default_domain_length")]
    pub domain_length: f64,
}

fn default_pressure() -> f64 { 101325.0 }
fn default_domain_length() -> f64 { 0.02 }

#[derive(Debug, Deserialize)]
pub struct GridConfig {
    #[serde(default = "default_initial_points")]
    pub initial_points: usize,
    #[serde(default = "default_max_points")]
    pub max_points: usize,
    #[serde(default = "default_grad")]
    pub grad: f64,
    #[serde(default = "default_curv")]
    pub curv: f64,
}

fn default_initial_points() -> usize { 20 }
fn default_max_points() -> usize { 500 }
fn default_grad() -> f64 { 0.05 }
fn default_curv() -> f64 { 0.10 }

#[derive(Debug, Deserialize)]
pub struct SolverConfig {
    #[serde(default = "default_atol")]
    pub atol: f64,
    #[serde(default = "default_rtol")]
    pub rtol: f64,
    #[serde(default = "default_max_newton")]
    pub max_newton_iter: usize,
    #[serde(default = "default_time_steps")]
    pub time_steps: usize,
    #[serde(default = "default_dt")]
    pub dt_initial: f64,
    pub su_initial_guess: Option<f64>,
    /// Optional path to a CSV initial profile (e.g. from Cantera).
    /// Overrides the sigmoid initial guess when provided.
    pub initial_profile: Option<String>,
}

fn default_atol() -> f64 { 1e-9 }
fn default_rtol() -> f64 { 1e-6 }
fn default_max_newton() -> usize { 50 }
fn default_time_steps() -> usize { 100 }
fn default_dt() -> f64 { 1e-7 }

#[derive(Debug, Deserialize)]
pub struct OutputConfig {
    #[serde(default = "default_output_file")]
    pub file: String,
}

fn default_output_file() -> String { "flame_solution.csv".to_string() }

impl FlameConfig {
    pub fn from_file(path: &str) -> Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let config: FlameConfig = toml::from_str(&content)?;
        config.validate()?;
        Ok(config)
    }

    fn validate(&self) -> Result<()> {
        anyhow::ensure!(self.flame.pressure > 0.0,
            "pressure must be positive, got {}", self.flame.pressure);
        anyhow::ensure!(self.flame.t_unburned > 0.0,
            "t_unburned must be positive, got {}", self.flame.t_unburned);
        anyhow::ensure!(self.flame.domain_length > 0.0,
            "domain_length must be positive, got {}", self.flame.domain_length);

        // Equivalence-ratio mode: all three fields required together.
        let has_fuel = self.flame.fuel.is_some();
        let has_ox   = self.flame.oxidizer.is_some();
        let has_phi  = self.flame.equivalence_ratio.is_some();
        let has_comp = self.flame.composition.is_some();

        anyhow::ensure!(
            has_comp ^ (has_fuel || has_ox || has_phi),
            "specify either [flame.composition] or all three of fuel/oxidizer/equivalence_ratio"
        );
        if !has_comp {
            anyhow::ensure!(has_fuel && has_ox && has_phi,
                "equivalence-ratio mode requires fuel, oxidizer, and equivalence_ratio");
            let phi = self.flame.equivalence_ratio.unwrap();
            anyhow::ensure!(phi > 0.0, "equivalence_ratio must be positive, got {phi}");
        }

        anyhow::ensure!(self.grid.initial_points >= 2,
            "grid.initial_points must be >= 2, got {}", self.grid.initial_points);
        anyhow::ensure!(self.grid.max_points >= self.grid.initial_points,
            "grid.max_points ({}) must be >= initial_points ({})",
            self.grid.max_points, self.grid.initial_points);
        anyhow::ensure!(self.grid.grad > 0.0 && self.grid.grad <= 1.0,
            "grid.grad must be in (0, 1], got {}", self.grid.grad);
        anyhow::ensure!(self.grid.curv > 0.0 && self.grid.curv <= 1.0,
            "grid.curv must be in (0, 1], got {}", self.grid.curv);

        anyhow::ensure!(self.solver.atol > 0.0,
            "solver.atol must be positive, got {}", self.solver.atol);
        anyhow::ensure!(self.solver.rtol > 0.0,
            "solver.rtol must be positive, got {}", self.solver.rtol);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_toml(extra_flame: &str) -> String {
        format!(r#"
[mechanism]
file = "data/h2o2.yaml"

[flame]
t_unburned = 300.0
{extra_flame}

[grid]

[solver]

[output]
"#)
    }

    fn parse(toml: &str) -> anyhow::Result<FlameConfig> {
        let cfg: FlameConfig = toml::from_str(toml)?;
        cfg.validate()?;
        Ok(cfg)
    }

    #[test]
    fn test_valid_composition_mode() {
        let toml = base_toml(r#"composition = { H2 = 0.3, O2 = 0.15, N2 = 0.55 }"#);
        assert!(parse(&toml).is_ok(), "valid composition mode should parse");
    }

    #[test]
    fn test_valid_equivalence_ratio_mode() {
        let toml = base_toml(r#"
equivalence_ratio = 1.0
[flame.fuel]
H2 = 1.0
[flame.oxidizer]
O2 = 0.21
N2 = 0.79
"#);
        assert!(parse(&toml).is_ok(), "valid φ mode should parse");
    }

    #[test]
    fn test_negative_pressure_rejected() {
        let toml = base_toml(r#"pressure = -1.0
composition = { N2 = 1.0 }"#);
        let err = parse(&toml).unwrap_err().to_string();
        assert!(err.contains("pressure"), "error should mention pressure: {err}");
    }

    #[test]
    fn test_zero_equivalence_ratio_rejected() {
        let toml = base_toml(r#"
equivalence_ratio = 0.0
[flame.fuel]
H2 = 1.0
[flame.oxidizer]
O2 = 0.21
N2 = 0.79
"#);
        let err = parse(&toml).unwrap_err().to_string();
        assert!(err.contains("equivalence_ratio"), "error should mention phi: {err}");
    }

    #[test]
    fn test_composition_and_phi_together_rejected() {
        let toml = base_toml(r#"
equivalence_ratio = 1.0
composition = { N2 = 1.0 }
[flame.fuel]
H2 = 1.0
[flame.oxidizer]
O2 = 0.21
N2 = 0.79
"#);
        assert!(parse(&toml).is_err(), "both modes at once should be rejected");
    }

    #[test]
    fn test_negative_domain_length_rejected() {
        let toml = base_toml(r#"domain_length = -0.01
composition = { N2 = 1.0 }"#);
        let err = parse(&toml).unwrap_err().to_string();
        assert!(err.contains("domain_length"), "error should mention domain_length: {err}");
    }
}
