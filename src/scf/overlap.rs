/*
Chance Brandt, The Ohio State Univeristy 2026

These are the gaussian integral solvers used to calculate the overlap integrals.

T. Helgaker, P. Taylor, 1995 "Gaussian Basis Sets and Molecular Integrals"
    - Hermite-recursive overlap of gaussians

"Fundamentals of Molecular Integrals Evaluation" by Justin T. Fermann and Edward F. Valeev 
    - Analytical solution to gaussian self overlap.
*/

use std::f64::consts::PI;
use factorial::{DoubleFactorial, Factorial};
use faer::traits::math_utils::sqrt;
use faer::Mat;
use crate::basis::{BasisSet, ContractedFunction};
use crate::molecule::Geometry;

#[inline]
pub fn hermite_coefficients(anga: &usize, expa: &f64, angb: &usize, expb: &f64, dist: &f64, nodes: &i8) -> f64 {
    /* Gets the expansion coefficients for Hermite polynomials recursively
    From Joshua Goings: https://joshuagoings.com/2017/04/28/integrals/
    
    anga/angb - orital angular momentum number of Gassian 'a' and 'b'
    expa/expb - orbital expontent on Gaussian 'a' and 'b'
    dist - distance between 
    
    TODO - make this not recursive so you can inline this shit
    */
    let p = expa + expb;
    let q = (expa * expb) / (p);

    if (*nodes < 0 || *nodes > (*anga + *angb) as i8) {
        // out of bounds
        return 0.0;
    } else if (*anga == 0 && *angb == 0 && *nodes == 0) {
        // base case
        return 1.0;
    } else if *angb == 0 {
        // decrement index a
        return (1.0/(2.0*p)) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &(nodes - 1)) - 
                ((q * dist) / expa) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &nodes) + 
                ((nodes + 1) as f64) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &(nodes + 1));

        //return (1.0/(2.0*p)) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &(nodes - 1));
        //return ((q * dist) / expa) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &nodes);
        //return ((nodes + 1) as f64) * hermite_coefficients(&(anga - 1), &expa, &angb, &expb, &dist, &(nodes + 1));
    } else {
        // decrement index b
        return (1.0/(2.0*p)) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &(nodes - 1)) + 
                ((q * dist) / expb) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &nodes) + 
                ((nodes + 1) as f64) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &(nodes + 1));
        
        //return (1.0/(2.0*p)) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &(nodes - 1));
        //return ((q * dist) / expb) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &nodes);
        //return((nodes + 1) as f64) * hermite_coefficients(&anga, &expa, &(angb - 1), &expb, &dist, &(nodes + 1));
    }
}

pub fn hermite_overlap(expa: &f64, angxyza: &[usize; 3], loca: &[f64; 3],
                       expb: &f64, angxyzb: &[usize; 3], locb: &[f64; 3]) -> f64{
    /* Gets the overlap integral between two gaussians using Hermite polynomials
    From Joshua Goings: https://joshuagoings.com/2017/04/28/integrals/
    */
    let S1 = hermite_coefficients(&angxyza[0], &expa, &angxyzb[0], &expb, &(loca[0] - locb[0]), &0); // X
    let S2 = hermite_coefficients(&angxyza[1], &expa, &angxyzb[1], &expb, &(loca[1] - locb[1]), &0); // Y
    let S3 = hermite_coefficients(&angxyza[2], &expa, &angxyzb[2], &expb, &(loca[2] - locb[2]), &0); // Z
    
    let p = expa + expb;
    let q = (expa * expb) / (p);
    let r2 = (loca[0] - locb[0]).powi(2) +
             (loca[1] - locb[1]).powi(2) + 
             (loca[2] - locb[2]).powi(2);

    return sqrt(&(PI/(p))).powi(3) * (-1.0 * q * r2).exp() * S1 * S2 * S3; 
}

pub fn hermite_contracted_overlap(coeffa: &[f64], expa: &[f64], angxyza: &[usize; 3], loca: &[f64; 3], prim_numa: &usize,
                                  coeffb: &[f64], expb: &[f64], angxyzb: &[usize; 3], locb: &[f64; 3], prim_numb: &usize) -> f64 {
    /* Gets the overlap integral between two congracted basis functions using Hermite
    polynomials.
    From Joshua Goings: https://joshuagoings.com/2017/04/28/integrals/
    
    Gets the overlap between two contracted gaussians
    */
    let mut S = 0.0;

    for i in 0..*prim_numa {
        for j in  0..*prim_numb {
            S += coeffa[i] * coeffb[j] * hermite_overlap(
                &expa[i], &angxyza, &loca, &expb[j], &angxyzb, &locb 
            );
        }
    }
    
    return S;
}

pub fn get_cartesian_overlap(basis: &BasisSet, mol: &Geometry) -> Mat::<f64> {
    /*
    Builds and returns a square matrix with the overlap between all contracted gaussians
    */    

    let mut overlap = Mat::<f64>::zeros(basis.total_aos, basis.total_aos);
    // This currently doesn't make use of the symmetric nature of the overlap matrix...
    // TODO
    let mut i = 0;
    for shelli in basis.shells.iter() {
        let loci = mol.coords[shelli.ele_offset];
        for funci in shelli.functions.iter() {
            let angxyzi = [funci.lx, funci.ly, funci.lz];
            let coeffi = &basis.prim_coeffs[funci.coeff_offset..(funci.coeff_offset + shelli.prim_num)];
            let expi = &basis.prim_exp[funci.exp_offset..(funci.exp_offset + shelli.prim_num)];
            let mut j = 0;
            for shellj in basis.shells.iter() {
                let locj = mol.coords[shellj.ele_offset];
                for funcj in shellj.functions.iter() {
                    let angxyzj = [funcj.lx, funcj.ly, funcj.lz];
                    let coeffj = &basis.prim_coeffs[funcj.coeff_offset..(funcj.coeff_offset + shellj.prim_num)];
                    let expj = &basis.prim_exp[funcj.exp_offset..(funcj.exp_offset + shellj.prim_num)];
                    let cur_overlap = hermite_contracted_overlap(
                        &coeffi, &expi, &angxyzi, &loci, &shelli.prim_num,
                        &coeffj, &expj, &angxyzj, &locj, &shellj.prim_num
                    );
                    overlap[(i, j)] += cur_overlap;
                    j += 1;
                }
            }
            i += 1;
        }
    }

    return overlap;
}

#[inline]
pub fn one_center_one_gaussian_integral(lx: usize, ly: usize, lz: usize, exp: f64, coeff: f64) -> f64 {
    /*  From "Fundamentals of Molecular Integrals Evaluation" by
            Justin T. Fermann and Edward F. Valeev 
        Analytical solution to gaussian self overlap.
    */
    let a = if (lx > 0) {2 * lx - 1} else {1};
    let b = if (ly > 0) {2 * ly - 1} else {1};
    let c = if (lz > 0) {2 * lz - 1} else {1};

    let L = lx + ly + lz;

    let integral = ((a.double_factorial() * b.double_factorial() * c.double_factorial()) as f64) * 
                    sqrt(&PI).powi(3) / 
                   ( (sqrt(&(2.0)) as f64).powi((L) as i32) * sqrt(&(2.0 * exp)).powi(3));
    //return integral * coeff * coeff;
    return integral;
}

#[inline]
pub fn one_center_two_gaussian_integral(lx: usize, ly: usize, lz: usize,
                                        exp1: f64, coeff1: f64,
                                        exp2: f64, coeff2: f64) -> f64
{
    /*  From "Fundamentals of Molecular Integrals Evaluation" by
            Justin T. Fermann and Edward F. Valeev 
        Analytical solution to gaussian overlap of two functions centered at the same point.
    */
    let a = if (lx > 0) {2 * lx - 1} else {1};
    let b = if (ly > 0) {2 * ly - 1} else {1};
    let c = if (lz > 0) {2 * lz - 1} else {1};

    let L = lx + ly + lz;

    let prefactor =  (sqrt(&PI.powi(3)) * 
        (a.double_factorial() * 
        b.double_factorial() * 
        c.double_factorial()) as f64) / 
        ((2 as i32).pow((L) as u32) as f64);

    return prefactor * ((coeff1 * coeff2) / ( sqrt(&(exp1 + exp2)).powi((2 * L + 3) as i32)));
}