/*
Chance Brandt, The Ohio State 2026

This is a toy quantum chemistry solver for learning, demonstrations,
and exploration features. It will be implemented in rust as much as
possible.

This is the main file used to orchestrate the rest of the program. It
starts with reading an input file and geometry file, and executing 
based on those instructions.

TODOs:

    - Print levels for handing different checks/debugging
        - Should come up with a scheme for what infor should be at each print level

    - Have defaults for all options? Could be tied in with the input variable obj
*/

mod file_input;
mod cli;
mod molecule;
mod scf;
mod print;
mod basis;
mod parse_json;
mod util;
mod testing;
use faer::Mat;

fn main() {

    let line_arg = cli::get_arguments().unwrap();
    let initial_data = file_input::read_input(line_arg.in_file.to_str()
                                                .expect("Cannot read input file"))
                                                .expect("Error in parsing geometry or input");

    // Should consider making this a callable object for easier handling throughout
    let section_info = initial_data.0;
    let mol = initial_data.1;
    
    mol.print();

    let mol_basis: basis::BasisSet = basis::load_basis(&mol, "STO3G.json");

    let overlap: Mat::<f64> = scf::overlap::get_cartesian_overlap(&mol_basis, &mol);
    print::mat_print::print_2D_mat(&overlap);

    println!("");
    let kinetic_mat: Mat::<f64> = scf::kinetic::get_kinetic_matrix(&mol_basis, &mol);
    print::mat_print::print_2D_mat(&kinetic_mat);

    println!("");
    let nuclear_mat: Mat::<f64> = scf::nuclear::get_nuclear_one_electron_matrix(&mol_basis, &mol);
    print::mat_print::print_2D_mat(&nuclear_mat);
}

