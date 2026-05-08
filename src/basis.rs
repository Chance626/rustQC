/*
Chance Brandt, The Ohio State University 2026

These are the structures and functions related to gettings and setting basis
functions from the ../basis/ objects and then structuring them into rust usable
objects.
*/

#[inline]
pub fn index_of(i: usize, j: usize) -> usize {
    /*
    Get's lower triangle matrix index of the ith and jth element
     */
    let (a, b) = if i > j {(i, j)} else {(j, i)};
    ( a * (a + 1) / 2) + b
}

pub enum BasisType {
    Cartesian,
    Spherical
}

#[repr(C)]
pub struct BasisSet {
    // Information about angular momement and contractions stored here
    pub shells: Vec<ContractedShell>,

    pub prim_coeffs: Vec<f64>,
    pub prim_exp: Vec<f64>
}

#[repr(C)]
pub struct ContractedShell {
    pub gauss_type: BasisType, 
    pub l: usize, // Angular Momentum

    // How to get the correct prim numbers from BasisSet
    pub prim_num: usize, // Number of primitives in the contraction
    pub coeff_offset: usize, // displacement from the start of prim_coeffs
    pub exp_offset: usize, // displacement from the start of prim_exp

    // How to get the correct coords from a Geometry
    pub ele_offset: usize, // displacement from the start of Geometry.eles

    // easy indexing to AO representation
    pub ao_offset: usize
}

impl BasisSet{
    pub fn normalize(&mut self) {
        // TODO: Add normalization
    }

    pub fn print() {
        // TODO: Add printing
    }
}

impl ContractedShell {
    #[inline]
    pub fn get_num_ao(&self) -> usize {
        match self.gauss_type {
            BasisType::Cartesian =>
                (self.l + 1) * (self.l + 2) / 2,
            
            BasisType::Spherical =>
                2 * self.l + 1,
        }
    }
}