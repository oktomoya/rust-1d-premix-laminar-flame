use crate::chemistry::species::Species;

#[derive(Debug, Clone)]
pub struct Arrhenius {
    pub a: f64,
    pub b: f64,
    pub ea: f64, // J/mol
}

#[derive(Debug, Clone)]
pub struct TroeParams {
    pub a: f64,
    pub t3: f64,
    pub t1: f64,
    pub t2: Option<f64>,
}

#[derive(Debug, Clone)]
pub enum RateType {
    Arrhenius { a: f64, b: f64, ea: f64 },
    Falloff {
        high: Arrhenius,
        low: Arrhenius,
        troe: Option<TroeParams>,
    },
}

#[derive(Debug, Clone)]
pub struct ThirdBodySpec {
    /// Per-species enhancement factors (species_index, efficiency)
    pub efficiencies: Vec<(usize, f64)>,
}

#[derive(Debug, Clone)]
pub struct Reaction {
    /// (species_index, stoichiometric coefficient)
    pub reactants: Vec<(usize, f64)>,
    pub products: Vec<(usize, f64)>,
    pub rate: RateType,
    pub reversible: bool,
    pub third_body: Option<ThirdBodySpec>,
}

/// A complete chemical mechanism: species list + reaction list.
#[derive(Debug, Clone)]
pub struct Mechanism {
    pub species: Vec<Species>,
    pub reactions: Vec<Reaction>,
}

impl Mechanism {
    pub fn n_species(&self) -> usize {
        self.species.len()
    }

    pub fn n_reactions(&self) -> usize {
        self.reactions.len()
    }

    pub fn species_index(&self, name: &str) -> Option<usize> {
        self.species.iter().position(|s| s.name == name)
    }
}
