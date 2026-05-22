/*
These are random functions I need to operators, and don't want to install a complete
library for some of these.
*/

#[inline]
pub fn index_of(i: usize, j: usize) -> usize {
    /*
    Get's lower triangle matrix index of the ith and jth element
     */
    let (a, b) = if i > j {(i, j)} else {(j, i)};
    ( a * (a + 1) / 2) + b
}