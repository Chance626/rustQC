/*
Chance Brandt, The Ohio State University 2026

These are methods and functions to get the nuclear Coulombic interactions
for a postive point charge and a negative gaussian distribution of charge

Based off of: https://joshuagoings.com/2017/04/28/integrals/
*/

use std::f64::consts::PI;
use crate::util::{boys, gaussian_product_center, l2_norm};
use crate::scf::overlap::hermite_coefficients;
use faer::Mat;
use crate::molecule::Geometry;
use crate::basis::BasisSet;

pub fn get_nuclear_one_electron_matrix(basis: &BasisSet, mol: &Geometry) -> Mat<f64> {
    /*
    This function takes a normalized basis set and a molecule object and returns 
    a matrix of the nuclear electron integrals, treating the electrons as a gaussian
    distrbution and the nuclei as point charges
    */

    let mut nuc_mat = Mat::<f64>::zeros(basis.total_aos, basis.total_aos);

    // Iterate through atomic function i,j then use the function pair to get a 
    // nuclear potential for each nuclei A
    let mut i = 0;
    for shelli in basis.shells.iter() {
        let loca = &mol.coords[shelli.ele_offset];
        let prim_numa = shelli.prim_num;

        for funci in shelli.functions.iter() {
            let coeffa = &basis.prim_coeffs[funci.coeff_offset..(funci.coeff_offset + prim_numa)];
            let expa = &basis.prim_exp[funci.exp_offset..(funci.exp_offset + prim_numa)];
            let angxyza = &[funci.lx, funci.ly, funci.lz];
            let mut j = 0;

            for shellj in basis.shells.iter() {
                let locb = &mol.coords[shellj.ele_offset];
                let prim_numb = shellj.prim_num;

                for funcj in shellj.functions.iter() {
                    let coeffb = &basis.prim_coeffs[funcj.coeff_offset..(funcj.coeff_offset + prim_numb)];
                    let expb = &basis.prim_exp[funcj.exp_offset..(funcj.exp_offset + prim_numb)];
                    let angxyzb = &[funcj.lx, funcj.ly, funcj.lz];

                    for atomi in 0..mol.natoms {
                        // Do the nuclear solver for each
                        let nuc_loc = &mol.coords[atomi];
                        let nuc_charge = &mol.eles[atomi];
                        
                        let coulomb_V = nuclear_contracted_gaussian(
                            coeffa, expa, angxyza, loca, &prim_numa,
                            coeffb, expb, angxyzb, locb, &prim_numb,
                            nuc_loc, nuc_charge
                        );

                        nuc_mat[(i, j)] += coulomb_V;
                    }
                    j += 1;
                }
            }

            i += 1;
        }
    }

    let nuc_mat = nuc_mat;
    return nuc_mat;
}

pub fn nuclear_contracted_gaussian(
    coeffa: &[f64], expa: &[f64], angxyza: &[usize; 3], loca: &[f64; 3], prim_numa: &usize,
    coeffb: &[f64], expb: &[f64], angxyzb: &[usize; 3], locb: &[f64; 3], prim_numb: &usize,
    nuc_loc: &[f64; 3], nuc_charge: &u8
) -> f64 {
    /*
    Runs over all primitives within a contracted function and returns the sum of the 
    nuclear-coulomb attraction integrals
     */
    let mut V_contracted = 0.0;
    for i in 0..*prim_numa {
        for j in 0..*prim_numb {

            V_contracted += nuclear_gaussian(
                &coeffa[i], &expa[i], angxyza, loca,
                &coeffb[i], &expb[i], angxyzb, locb,
                nuc_loc, nuc_charge
            );
        }
    }
    return V_contracted;
}

pub fn nuclear_gaussian(
    coeffa: &f64, expa: &f64, angxyza: &[usize; 3], loca: &[f64; 3],
    coeffb: &f64, expb: &f64, angxyzb: &[usize; 3], locb: &[f64; 3],
    nuc_loc: &[f64; 3], nuc_charge: &u8
) -> f64 {
    /*
    Gets the nuclear-gaussian charge energy for a gaussian pair of orbitals and a
    point charge
     */
    let max_ord_x = (angxyza[0] + angxyzb[0]) as i32;
    let max_ord_y = (angxyza[1] + angxyzb[1]) as i32;
    let max_ord_z = (angxyza[2] + angxyzb[2]) as i32;

    let exp_sum = expa + expb;
    let gauss_center = gaussian_product_center(expa, loca, expb, locb);
    let gauss_point_dist = l2_norm(&[
        gauss_center[0] - nuc_loc[0],
        gauss_center[1] - nuc_loc[1],
        gauss_center[2] - nuc_loc[2]
    ]);

    let mut V = 0.0;
    for i in 0..(max_ord_x + 1) {
        for j in 0..(max_ord_y + 1) {
            for k in 0..(max_ord_z + 1) {
                V += (*nuc_charge as f64) *
                hermite_coefficients(&angxyza[0], expa, &angxyzb[0], expb, &(loca[0] - locb[0]), &(i as i8)) *
                hermite_coefficients(&angxyza[1], expa, &angxyzb[1], expb, &(loca[1] - locb[1]), &(j as i8)) *
                hermite_coefficients(&angxyza[2], expa, &angxyzb[2], expb, &(loca[2] - locb[2]), &(k as i8)) *
                coulomb_aux_hermite_integral(&i, &j, &k, &0, &exp_sum, &(gauss_center[0] - nuc_loc[0]), &(gauss_center[1] - nuc_loc[1]), &(gauss_center[2] - nuc_loc[2]), &gauss_point_dist);
            }
        }
    }
    V *= 2.0 * PI / exp_sum;
    return V;
}

pub fn coulomb_aux_hermite_integral(ord_x: &i32, ord_y: &i32, ord_z: &i32, n: &i32, exp_sum: &f64,
                                    dist_x: &f64, dist_y: &f64, dist_z: &f64, dist: &f64) -> f64 {
    /*
    Evaluates the Coulomb auxiliary Hermite integrals recursively.
    */
    let T: f64 = exp_sum * dist * dist;
    let mut val: f64 = 0.0;
    
    if (*ord_x == 0 && *ord_y == 0 && *ord_z == 0) {
        // Base Case
        val += (-2.0 * exp_sum).powi(*n) * boys(*n as f64, T);
    } else if (*ord_x == *ord_y && *ord_x == 0) {
        // Deal with ord_z
        if (*ord_z > 1) {
            val += (ord_z - 1) as f64 * coulomb_aux_hermite_integral(
                ord_x, ord_y, &(ord_z - 2), &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
            );
        }
        val += dist_z * coulomb_aux_hermite_integral(
            ord_x, ord_y, &(ord_z - 1), &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
        );
    } else if (*ord_x == 0) {
        // Deal with ord_y
        if (*ord_y > 1) {
            val += (ord_y - 1) as f64 * coulomb_aux_hermite_integral(
                ord_x, &(ord_y -2), ord_z, &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
            );
        }
        val += dist_y * coulomb_aux_hermite_integral(
            ord_x, &(ord_y - 1), ord_z, &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
        );

    } else {
        // Deal with ord_x
        if (*ord_x > 1) {
            val += (ord_x - 1) as f64 * coulomb_aux_hermite_integral(
                &(ord_x - 2), ord_y, ord_z, &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
            );
            val += dist_x * coulomb_aux_hermite_integral(
                &(ord_x - 1), ord_y, ord_z, &(n + 1), exp_sum, dist_x, dist_y, dist_z, dist
            );
        }
    }

    return val;
}