/*
These are random functions I need to operators, and don't want to install a complete
library for some of these.
*/

use factorial::Factorial;
use numra_special::hypergeometric::hyp1f1;

#[inline]
pub fn index_of(i: usize, j: usize) -> usize {
    /*
    Get's lower triangle matrix index of the ith and jth element
     */
    let (a, b) = if i > j {(i, j)} else {(j, i)};
    ( a * (a + 1) / 2) + b
}

#[inline]
pub fn binom(n: u32, k: u32) -> u32 {
    return n.factorial() / (k.factorial() * (n - k).factorial());
}

#[inline]
pub fn boys(n: f64, T: f64) -> f64 {
    /*
    Wrapper for the confluent hypergeometric function to account for the 
    case of:
    F_n(T) = int_{0}^{1} ( exp(-Tx^2)x^(2n)dx )
    Called the boy's function
     */
    return hyp1f1(n+0.5, n+1.5, -T);
}

#[inline]
pub fn gaussian_product_center(expA: &f64, locA: &[f64; 3], expB: &f64, locB: &[f64; 3]) -> [f64; 3] {
    let mut prod_center = [0.0; 3];
    let denom = (expA + expB);
    for i in 0..3 {
        prod_center[i] = ((expA * locA[i]) + (expB * locB[i])) / denom;
    }
    let prod_center = prod_center;
    return prod_center;
}

#[inline]
pub fn l2_norm(curvec: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..curvec.len() {
        sum += curvec[i] * curvec[i];
    }
    return sum.sqrt();
}