/*
Chance Brandt, The Ohio State University 2026

This file handles input file parsing and stores the user variables
and default values for the rest of the run. Largely controls the
initialization of the run.
*/

pub fn read_input() {
    // Reads and stores the input setting portion of the input settings
    // section and stores them for future use in the program.

    println!("Starting in read_input().");
    read_geometry();
    println!("Ending in read_input().");

}

pub fn read_geometry() {
    // Reads and stores the geometry portion of the input file or 
    // the specified geometry file. Currently can accept:
    //  .xyz

    println!("Hello from read_geometry().");

}