/*
Chance Brandt, The Ohio State University 2026

These are the structures and functions related to gettings and setting basis
functions from the ../basis/ objects and then structuring them into rust usable
objects.
*/

use std::{f64::consts::PI, fs};
use crate::{molecule, parse_json};
use std::cmp::{PartialEq, Eq};
use faer::{self, traits::math_utils::sqrt};
use crate::scf::integral_solver::single_integral;

#[derive(PartialEq, Eq, Debug)]
pub enum BasisType {
    Cartesian,
    Spherical
    // TODO: implement switching between string and BasisType enum, jank if/else below
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
        // This normalizes all the coefficients

    }

    pub fn get_prim_norms(&self) -> Vec<f64> {

        let mut sum_prim = 0;

        for shell in self.shells.iter() {
            sum_prim += shell.prim_num;
        }

        let mut norms: Vec<f64> = vec![0.0; sum_prim];
        let mut count = 0;
        for shell in self.shells.iter() {
            for i in 0..shell.prim_num {
                let cur_coef = self.prim_coeffs[shell.coeff_offset + i];
                let cur_exp = self.prim_exp[shell.exp_offset + i];
                let cur_overlap = single_integral(shell.l, cur_exp, cur_coef);
                // ensures that the norms and exps are indexed the same
                norms[shell.exp_offset + i] = (1.0 / cur_overlap);
                count += 1;
            }
        }

        let norms = norms;
        return norms;
    }

    pub fn print(&self) {
        let header = "> Basis Sets <";
        println!("{:=^48}\n", header);
        for shell in self.shells.iter() {
            println!("{:7}{:?} Shell for atom: {}, l = {}",
                "",
                shell.gauss_type,
                shell.ele_offset + 1,
                shell.l);
            println!("{:9}Coefficients{:8}Exponents",
            "",
            "");
            for i in 0..shell.prim_num {
                println!("{:11}{:<12.9}{:8}{:<12.9}",
                    "",
                    self.prim_coeffs[i + shell.coeff_offset],
                    "",
                    self.prim_exp[i + shell.exp_offset]
                );
            }
        }
        println!("{:=^48}\n", "");
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

pub fn load_basis (geom: &molecule::Geometry, basis_name: &str) -> BasisSet {
    let base_path = std::env!("CARGO_MANIFEST_DIR");
    let basis_path = &format!("{base_path}/basis/{basis_name}");
    let contents = &fs::read_to_string(basis_path).unwrap();

    let ele_basis_sets: parse_json::json_basis = serde_json::from_str(contents).unwrap();

    let mol_basis: BasisSet = build_mol_basis(&geom, &ele_basis_sets, BasisType::Cartesian);

    return mol_basis;
}

pub fn build_mol_basis (geom: &molecule::Geometry, ele_basis_sets: &parse_json::json_basis, mode: BasisType) -> BasisSet {


    let func_type: &str = if mode == BasisType::Spherical {"gto_spherical"} else {"gto"};
    let mut coef_offset = 0 ;
    let mut exp_offset = 0 ;
    let mut ao_offset = 0 ;
    let mut mol_basis = BasisSet { shells: Vec::new(), prim_coeffs: Vec::new(), prim_exp: Vec::new() };
    
    for i in 0..geom.natoms {
        let cur_ele: u8 = geom.eles[i];
        let cur_shells: &Vec<parse_json::electron_shell> = &ele_basis_sets.elements[&(cur_ele as usize)].electron_shells;
        for shell in cur_shells.iter() {
            if shell.function_type == func_type {
                mol_basis.prim_exp.extend(shell.exponents.iter().copied());
                for ang in shell.angular_momentum.iter() {
                    let cur_basis_type: BasisType = if func_type == "gto_spherical" {BasisType::Spherical} else {BasisType::Cartesian};
                    mol_basis.prim_coeffs.extend(shell.coefficients[*ang].iter().copied());
                    let cur_contract_shell: ContractedShell = 
                        ContractedShell {
                            gauss_type: cur_basis_type,
                            l: *ang,
                            prim_num: shell.exponents.len(),
                            coeff_offset: coef_offset,
                            exp_offset: exp_offset,
                            ele_offset: i,
                            ao_offset: ao_offset
                        };
                        coef_offset += shell.coefficients[*ang].len();
                        ao_offset += cur_contract_shell.get_num_ao();
                        mol_basis.shells.push(cur_contract_shell);
                }
                exp_offset += shell.exponents.len();
            }
        }
    }

    let mol_basis = mol_basis;
    return mol_basis;
}