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
pub fn self_integral(lx: usize, ly: usize, lz: usize, exp: f64, coeff: f64) -> f64 {
    /* From "Fundamentals of Molecular Integrals Evaluation" by
            Justin T. Fermann and Edward F. Valeev */
    let a = if (lx > 0) {2 * lx - 1} else {1};
    let b = if (lx > 0) {2 * lx - 1} else {1};
    let c = if (lx > 0) {2 * lx - 1} else {1};

    let integral = (((a.double_factorial() * b.double_factorial() * c.double_factorial()) as f64) * 
                    sqrt(&PI)).powi(3) / 
                   ((4.0 * exp).powi((lx + ly + lz) as i32) * sqrt(&(2.0 * exp)).powi(3));
    return integral * coeff * coeff;
}
