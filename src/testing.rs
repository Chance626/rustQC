/*
Chance Brandt, The Ohio State University 2026

Random functions I use for testing various parts of code, usually run from main.
*/

use crate::basis::BasisSet;
use crate::molecule;
use crate::scf::overlap;
use faer::Mat;

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
            println!("Hermite primitive overlap = {}", overlap::hermite_contracted_overlap(
                &cur_coeff, &cur_exp, &angxyz, &cur_loc, &shell.prim_num,
                &cur_coeff, &cur_exp, &angxyz, &cur_loc, &shell.prim_num));
        }
    }

}
