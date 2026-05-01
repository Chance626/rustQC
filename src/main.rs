/*
Chance Brandt, The Ohio State 2026

This is a toy quantum chemistry solver for learning, demonstrations,
and exploration features. It will be implemented in rust as much as
possible.

This is the main file used to orchestrate the rest of the program. It
starts with reading an input file and geometry file, and executing 
based on those instructions.
*/

mod file_input;
mod cli;

fn main() {

    cli::get_arguments();
    file_input::read_input("/users/PAS0291/cbrandt/tests/rust/rustQC/input.inp".to_string());

}