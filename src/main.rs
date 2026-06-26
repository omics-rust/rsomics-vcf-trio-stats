mod cli;

use std::process::ExitCode;

use clap::Parser;

use cli::Cli;

fn main() -> ExitCode {
    let args = Cli::parse();
    let json = args.common.json;
    match args.execute() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            if json {
                let env = serde_json::json!({
                    "schema_version": rsomics_common::SCHEMA_VERSION,
                    "tool": cli::META.name,
                    "tool_version": cli::META.version,
                    "status": "error",
                    "error": { "message": e.to_string() },
                });
                eprintln!("{}", serde_json::to_string(&env).unwrap_or_default());
            } else {
                eprintln!("error: {e}");
            }
            ExitCode::FAILURE
        }
    }
}
