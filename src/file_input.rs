/*
Chance Brandt, The Ohio State University 2026

This file handles input file parsing and stores the user variables
and default values for the rest of the run. Largely controls the
initialization of the run.
*/

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

pub fn read_input(filepath: String) -> i8 {
    // Reads and stores the input setting portion of the input settings
    // section and stores them for future use in the program.


    let section_info: HashMap<String, String> = read_sections(filepath);

    for (k, v) in section_info.iter() {
        println!("{} contains:\n{}", k, v)
    }

    read_geometry();

    return 0;
}

fn read_sections(filepath: String) -> HashMap<String, String> {
    /*
    Gets each section of the input file from the '%string' character to the next
    '%end', all information will be saved in a hashmap with the key of 'string'
    and the informations between the controls characters as the values.

    section information is returned with no blank lines or trailing/preceeding
    tabs/spaces

    comments are denoted by the '#' character
     */
    let mut sections = HashMap::new();
    
    // Error handling for the file open
    let file = match File::open(filepath.clone()) {
        Ok(f) => f,
        Err(e) => {
            eprint!("{e} for the file: {filepath}");
            return HashMap::new();
        }
    };
    let reader = BufReader::new(file);

    let mut section_contents = String::new();
    let mut section_header = String::new();
    let mut in_section: bool = false;

    for line in reader.lines() {
        let mut cur_line: String = line.unwrap().trim().to_string();
        if cur_line.starts_with("#") {
            continue;
        }
        
        if cur_line.starts_with("%") {
            cur_line = cur_line.trim_start_matches("%").to_string();
            if cur_line.to_lowercase() == "end" {
                in_section = false;
                sections.insert(section_header.clone(), section_contents.clone());
                
                section_header.clear();
                section_contents.clear();
            } else {
                in_section = true;
                section_header = cur_line;
            }
        } else if in_section && !cur_line.is_empty() {
            section_contents.push_str(&cur_line);
            section_contents.push('\n');
        }
    }

    return sections;
}

pub fn read_geometry() {
    /*
    Reads and stores the geometry portion of the input file or
    the specified geometry file. Currently can accept:
    .xyz
     */

    println!("Hello from read_geometry().");

}

