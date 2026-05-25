
/*
Chance Brandt, The Ohio State 2026

This is a toy quantum chemistry solver for learning, demonstrations,
and exploration features. It will be implemented in rust as much as
possible.

This is the main file used to orchestrate the rest of the program. It
starts with reading an input file and geometry file, and executing 
based on those instructions.
*/

mod file_input;
mod cli;
mod molecule;
mod scf;
mod print;
mod basis;
mod parse_json;
mod util;
// mod context;

fn main() {

    let line_arg = cli::get_arguments().unwrap();
    let initial_data = file_input::read_input(line_arg.in_file.to_str()
                                                .expect("Cannot read input file"))
                                                .expect("Error in parsing geometry or input");

    let section_info = initial_data.0;
    let mol = initial_data.1;
    
    mol.print();

    // need to perform operations on the basis before leaving it alone for the 
    // rest of the run, should abstract at some point
    let mut mol_basis: basis::BasisSet = basis::load_basis(&mol, "STO3G.json");
    mol_basis.print();
    mol_basis.normalize();
    //mol_basis.print();

    let mol_basis = mol_basis;
    // run the method of the input file
}