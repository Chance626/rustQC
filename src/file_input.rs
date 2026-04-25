/*
Chance Brandt, The Ohio State University 2026

This file handles input file parsing and stores the user variables
and default values for the rest of the run. Largely controls the
initialization of the run.
*/

use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn read_input(filepath: String) {
    // Reads and stores the input setting portion of the input settings
    // section and stores them for future use in the program.


    read_sections(filepath);

    read_geometry();

}

fn read_sections(filepath: String) {
    // Error handling for the file open
    let file = match File::open(filepath.clone()) {
        Ok(f) => f,
        Err(e) => {
            eprint!("{e} for the file: {filepath}");
            return;
        }
    };
    let reader = BufReader::new(file);

    //reader.lines() exists within 
    for cur_line in reader.lines() {
        // janked the error handling for once in the file
        println!("{}", cur_line.unwrap());
    }

}

pub fn read_geometry() {
    // Reads and stores the geometry portion of the input file or 
    // the specified geometry file. Currently can accept:
    //  .xyz

    println!("Hello from read_geometry().");

}

