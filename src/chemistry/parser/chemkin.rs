/// Parser for CHEMKIN-II format:
///   - chem.inp: element/species/reaction declarations
///   - therm.dat: NASA polynomial thermodynamic data
///   - tran.dat:  transport parameters
///
/// This is a skeletal implementation; full CHEMKIN parsing is extensive.
/// For now it provides the entry-point signatures and a minimal line tokenizer.

use anyhow::{bail, Result};
use crate::chemistry::mechanism::Mechanism;

pub struct ChemkinFiles {
    pub chem_inp: String,
    pub therm_dat: String,
    pub tran_dat: String,
}

/// Parse a CHEMKIN-II mechanism from text content.
pub fn parse(files: &ChemkinFiles) -> Result<Mechanism> {
    let _elements = parse_elements(&files.chem_inp)?;
    let species_names = parse_species_list(&files.chem_inp)?;
    let species = parse_thermo(&files.therm_dat, &species_names)?;
    let species = attach_transport(species, &files.tran_dat)?;
    let reactions = parse_reactions(&files.chem_inp, &species)?;

    Ok(Mechanism { species, reactions })
}

fn parse_elements(src: &str) -> Result<Vec<String>> {
    // TODO: extract ELEMENTS ... END block
    let _ = src;
    Ok(vec![])
}

fn parse_species_list(src: &str) -> Result<Vec<String>> {
    // TODO: extract SPECIES ... END block
    let _ = src;
    Ok(vec![])
}

fn parse_thermo(
    therm: &str,
    species_names: &[String],
) -> Result<Vec<crate::chemistry::species::Species>> {
    // TODO: parse NASA 7-coeff blocks
    // Format: species_name, date, elements, phase, T_low T_high T_mid
    //         5 coeffs (high range), 5 more, 4+2 last
    let _ = (therm, species_names);
    bail!("CHEMKIN thermo parser not yet implemented")
}

fn attach_transport(
    mut species: Vec<crate::chemistry::species::Species>,
    tran: &str,
) -> Result<Vec<crate::chemistry::species::Species>> {
    // TODO: parse tran.dat (one line per species)
    let _ = tran;
    Ok(species)
}

fn parse_reactions(
    src: &str,
    species: &[crate::chemistry::species::Species],
) -> Result<Vec<crate::chemistry::mechanism::Reaction>> {
    // TODO: parse REACTIONS ... END block
    //   - Arrhenius A, b, Ea on same line as stoichiometry
    //   - Optional LOW/, TROE/, PLOG/, DUPLICATE keywords
    let _ = (src, species);
    bail!("CHEMKIN reaction parser not yet implemented")
}
