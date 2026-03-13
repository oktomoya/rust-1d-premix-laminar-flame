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
        Ok(config)
    }
}
