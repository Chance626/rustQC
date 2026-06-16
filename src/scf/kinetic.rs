/*
Chance Brandt, The Ohio State University 2026

These are the methods and functions used to solve for the one-electron momentum
and kinect energy operators

T. Helgaker, P. Taylor, 1995 "Gaussian Basis Sets and Molecular Integrals"
    - Hermite-recursive overlap of kinetic matrices
*/

use faer::Mat;

use crate::scf::overlap;
use crate::basis::BasisSet;
use crate::molecule::Geometry;

pub fn get_kinetic_matrix(basis: &BasisSet, mol: &Geometry) -> Mat<f64> {
    /*
    Takes a normalized basis set and a molecule object and returns a matrix of 
    the kinectic energy values between the various orbitals.
    */
    let mut kin_mat = Mat::<f64>::zeros(basis.total_aos, basis.total_aos);
    
    let mut i = 0;
    for shelli in basis.shells.iter() {
        let loca = mol.coords[shelli.ele_offset];
        let prim_numa = shelli.prim_num;
        for funci in shelli.functions.iter() {
            let coeffa = &basis.prim_coeffs[funci.coeff_offset..(funci.coeff_offset + prim_numa)];
            let expa = &basis.prim_exp[funci.exp_offset..(funci.exp_offset + prim_numa)];
            let angxyza = [funci.lx, funci.ly, funci.lz];
            let mut j = 0;
            for shellj in basis.shells.iter() {
                let locb = mol.coords[shellj.ele_offset];
                let prim_numb = shellj.prim_num;
                for funcj in shellj.functions.iter() {
                    let coeffb = &basis.prim_coeffs[funcj.coeff_offset..(funcj.coeff_offset + prim_numb)];
                    let expb = &basis.prim_exp[funcj.exp_offset..(funcj.exp_offset + prim_numb)];
                    let angxyzb = [funcj.lx, funcj.ly, funcj.lz];
                    let cur_kin: f64 = contracted_kinetic(
                        coeffa, expa, &angxyza, &loca, &prim_numa,
                        coeffb, expb, &angxyzb, &locb, &prim_numb
                    );

                    kin_mat[(i, j)] += cur_kin;
                    j += 1;
                }
            }

            i += 1;
        }
    }

    return kin_mat;
}

pub fn contracted_kinetic(coeffa: &[f64], expa: &[f64], angxyza: &[usize; 3], loca: &[f64; 3], prim_numa: &usize,
               coeffb: &[f64], expb: &[f64], angxyzb: &[usize; 3], locb: &[f64; 3], prim_numb: &usize) -> f64 {
    /*
    Takes two contracted gaussians as parameters and returns the kinetic energy 
    integral between the two gaussians. (One electron, two center integral)

    T. Helgaker, P. Taylor, 1995 "Gaussian Basis Sets and Molecular Integrals"
    - Hermite-recursive overlap of kinetic matrices
    - section 11.2

    This sums the contracted gaussian components for a given pair of shells.
    */

    let mut kinetic_sum = 0.0;

    for i in 0..*prim_numa {
        for j in 0..*prim_numb {
            kinetic_sum += -0.5 * (coeffa[i] * coeffb[j] * laplacian_one_electron_two_center(
                &expa[i], angxyza, loca, &expb[j], angxyzb, locb
            ));
        }
    }

    return kinetic_sum;
}

pub fn laplacian_one_electron_two_center(expa: &f64, angxyza: &[usize; 3], loca: &[f64; 3],
                                         expb: &f64, angxyzb: &[usize; 3], locb: &[f64; 3]) -> f64 {
    /*
    From Taylor 1995, returns the D^2 terms of the kinetic energy term, which is 
    the second spatial derivative of the overlap (i.e. <Ga| d^2/dx^2 |Gb>), which
    are spatially separable.
    
    This term is expressed as a function of the angular momentum per cartesian 
    coordinates and the zero-th order overlap:
    
    D^2_i,j = j (j - 1) (S^0_i, j - 2) - 2b (2j + 1) (S^0_i,j) + 4b^2 (S^0_i,j + 2) 
    */
    
    let mut cart_summed_laplacian = 0.0;
    
    // summing over cartesian coordinates
    for i in 0..3 {
        // Decrement each angb by 2
        let mut dec_angxyzb = angxyzb.clone();
        dec_angxyzb[i] = dec_angxyzb[i].saturating_sub(2 as usize);
        let term1 = (angxyzb[i] as f64 * (angxyzb[i] as f64 - 1.0)) * overlap::hermite_overlap(expa, angxyza, loca, expb, &dec_angxyzb, locb) ;
        
        // No changes in angular momentum
        let term2 = 2.0 * expb * (2 * angxyzb[i] + 1) as f64 * overlap::hermite_overlap(expa, angxyza, loca, expb, angxyzb, locb);
        
        // Increment each angb by 2
        let mut inc_angxyzb = angxyzb.clone();
        inc_angxyzb[i] += 2;
        let term3 = 4.0 * (expb * expb) * overlap::hermite_overlap(expa, angxyza, loca, expb, &inc_angxyzb, locb);
        cart_summed_laplacian += term1 - term2 + term3;
    }

    return cart_summed_laplacian;
}