
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

use crate::{basis::BasisSet, scf::integral_solver};

mod file_input;
mod cli;
mod molecule;
mod scf;
mod print;
mod basis;
mod parse_json;
mod util;

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

    let mut mol_basis: basis::BasisSet = basis::load_basis(&mol, "STO3G.json");
    let mol_basis = mol_basis;


    let overlap: Mat::<f64> = integral_solver::get_cartesian_overlap(&mol_basis, &mol);
    print::mat_print::print_2D_mat(&overlap);
}

/*
// was testing analytical vs. hermite overlaps
pub fn testing_1(mol_basis: &BasisSet, mol: &molecule::Geometry, overlap: &Mat::<f64>) {

    mol_basis.print_contracted_self_overlap("Unnormalized Contracted Self_Overlap");

    for shell in mol_basis.shells.iter() {
        for func in shell.functions.iter() {
            let cur_norm = mol_basis.contract_norms[func.coeff_offset];
            let cur_norm_overlap = 1.0 / (cur_norm * cur_norm);
            let cur_exp = &mol_basis.prim_exp[func.exp_offset..(func.exp_offset + shell.prim_num)];
            let cur_coeff = &mol_basis.prim_coeffs[func.coeff_offset..(func.coeff_offset + shell.prim_num)];
            let cur_loc = mol.coords[shell.ele_offset];
            let angxyz = [func.lx, func.ly, func.lz];
            println!("atom: {}, lx: {}, ly: {}, lz: {}", 
                    shell.ele_offset, func.lx, func.ly, func.lz);
            //println!("Analytical norm overlap = {}", cur_norm_overlap);
            println!("Hermite primitive overlap = {}", integral_solver::hermite_contracted_overlap(
                &cur_coeff, &cur_exp, &angxyz, &cur_loc, &shell.prim_num,
                &cur_coeff, &cur_exp, &angxyz, &cur_loc, &shell.prim_num));
        }
    }

}

*/