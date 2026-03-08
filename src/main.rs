use anyhow::Result;
use clap::Parser;
use premix1d::io::input::FlameConfig;
use premix1d::flame::solver_driver::run_flame;

#[derive(Parser)]
#[command(name = "premix1d", about = "1D freely propagating premixed laminar flame solver")]
struct Args {
    /// Path to TOML input file
    #[arg(short, long, default_value = "flame.toml")]
    input: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let config = FlameConfig::from_file(&args.input)?;
    run_flame(&config)?;
    Ok(())
}
