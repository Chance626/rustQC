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
pub fn single_integral(ang: usize, exp: f64, coeff: f64) -> f64 {
    let a = if (ang > 0) {2 * ang - 1} else {1};
    let integral = (((a.double_factorial() as f64)) * sqrt(&PI)) / (  (4.0 * exp).powi(ang as i32) * sqrt(&(2.0 * exp)));
    return integral * coeff * coeff;
}
