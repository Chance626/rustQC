/*
Chance Brandt, The Ohio State Univeristy 2026

These are the gaussian integral solvers used to calculate the overlap integrals.

Let see how this goes lol
*/

use std::f64::consts::PI;
use factorial::DoubleFactorial;
use faer::traits::math_utils::sqrt;

#[inline]
pub fn E() {
    /* Gets the expansion coefficients for Hermite polynomials recursively*/
    
}

#[inline]
pub fn one_center_one_gaussian_integral(lx: usize, ly: usize, lz: usize, exp: f64, coeff: f64) -> f64 {
    /* From "Fundamentals of Molecular Integrals Evaluation" by
            Justin T. Fermann and Edward F. Valeev */
    let a = if (lx > 0) {2 * lx - 1} else {1};
    let b = if (ly > 0) {2 * ly - 1} else {1};
    let c = if (lz > 0) {2 * lz - 1} else {1};

    let integral = (((a.double_factorial() * b.double_factorial() * c.double_factorial()) as f64) * 
                    sqrt(&PI)).powi(3) / 
                   ((4.0 * exp).powi((lx + ly + lz) as i32) * sqrt(&(2.0 * exp)).powi(3));
    return integral * coeff * coeff;
}

#[inline]
pub fn one_center_two_gaussian_integral(lx: usize, ly: usize, lz: usize, 
                                        exp1: f64, coeff1: f64,
                                        exp2: f64, coeff2: f64) -> f64
{
    /* From "Fundamentals of Molecular Integrals Evaluation" by
            Justin T. Fermann and Edward F. Valeev 
            This doesn't apply the x, y, z prefactors so that the double
            factorials only have to be called once per contracted function*/
    let L = lx + ly + lz;
    return (coeff1 * coeff2) / ( sqrt(&(exp1 + exp2).powi((L + 3) as i32)));
}