/*
    Functions and structures for command line interfacing
*/

use std::env;
use std::path;

pub struct CliArgs{
    /* 
    Stores command line arguments for run execution
    */
    pub in_file: path::PathBuf,
    pub out_file: path::PathBuf,

}

pub fn get_arguments() -> Result<CliArgs, String> {
    let in_file = env::args().nth(1).expect("No in_file detected.");
    let out_file = env::args().nth(2).unwrap_or("output.out".to_string());
    let args = CliArgs{
        in_file: path::PathBuf::from(in_file),
        out_file: path::PathBuf::from(out_file)
    };

    Ok(args)
}