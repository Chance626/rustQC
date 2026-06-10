/*
These are random functions I need to operators, and don't want to install a complete
library for some of these.
*/

use factorial::Factorial;

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