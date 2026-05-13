/*
Chance Brandt, The Ohio State 2026

This is a wrapper function to be able to pass context information around a job run

Not currently implemented or structured properly
*/

mod basis;
mod molecule;

pub struct Context {
    // These are 
    pub geometry: &'a molecule::Geometry,
    pub basis: &'a basis::BasisSet,

}