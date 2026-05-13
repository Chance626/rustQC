/*
Chance Brandt, The Ohio State University 2026

This file handles input file parsing and stores the user variables
and default values for the rest of the run. Largely controls the
initialization of the run.
*/

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

use crate::molecule;

pub fn read_input(filepath: &str) -> Result<(HashMap<String, String>,
                                               molecule::Geometry), String> {
    /* 
    Reads and stores the input setting portion of the input settings
    section and stores them for future use in the program.
    */

    let section_info: HashMap<String, String> = read_sections(filepath);

    //for (k, v) in section_info.iter() {
    //    println!("{} contains:\n{}", k, v)
    //}

    let geom_name = String::from("geom");
    let mol: molecule::Geometry = read_geometry(section_info.get(&geom_name)
    .expect("Cannot read geom section."))?;

    Ok((section_info, mol))
}

fn read_sections(filepath: &str) -> HashMap<String, String> {
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
    let file = match File::open(filepath) {
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

fn read_geometry(geom: &str) -> Result<molecule::Geometry, String> {
    /*
    Reads and stores the geometry portion of the input file or
    the specified geometry file. Currently can accept:
    - geom input sections
     */

    let geom_lines: Vec<&str> = geom.lines().collect();
    let natoms: usize = geom_lines.len() - 1;
    
    let chrg_mult_line: Vec<&str> = geom_lines[0].split_whitespace().collect();
    let chrg: i8 = chrg_mult_line[0].parse()
    .map_err(|_| "Cannot read charge value")?;
    let mult: i8 = chrg_mult_line[1].parse()
    .map_err(|_| "Cannot read multiplicity value")?;
    let mut eles: Vec<u8> = vec![0; natoms];
    let mut coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0]; natoms];

    for i in 0..(geom_lines.len() - 1) {
        let cur_line: Vec<&str> = geom_lines[i + 1].split_whitespace().collect();

        let atom_num = molecule::ELEMENTS.get(cur_line[0])
        .ok_or_else(|| format!("Cannot identify element {}", cur_line[0]))?;
        
        eles[i] = *atom_num;

        for j in 0..3 {
            let cur_coords: f64 = cur_line[j + 1].parse()
            .map_err(|e| format!("Cannot interpret geometry coordinates {}", e))?;
            coords[i][j] = cur_coords;
        }
    }


    Ok(molecule::Geometry{
        eles: eles,
        coords: coords,
        natoms: natoms,
        chrg: chrg,
        mult: mult
    })
}