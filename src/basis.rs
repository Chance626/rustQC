/*
Chance Brandt, The Ohio State University 2026

These are the structures and functions related to gettings and setting basis
functions from the ../basis/ objects and then structuring them into rust usable
objects.
*/

use core::num;
use std::{f64::consts::PI, fs, primitive};
use crate::{molecule, parse_json};
use std::cmp::{PartialEq, Eq};
use faer::{self, traits::math_utils::sqrt};
use crate::scf::overlap::{
    one_center_one_gaussian_integral,
    one_center_two_gaussian_integral};
use factorial::DoubleFactorial;

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
    pub prim_exp: Vec<f64>,

    // These are the primitive and contractions normalization constants
    pub prim_norms: Vec<f64>,
    pub contract_norms: Vec<f64>,

    pub total_aos: usize
}

#[repr(C)]
pub struct ContractedShell {
    pub gauss_type: BasisType,
    pub l: usize, // Angular Momentum
    pub functions: Vec<ContractedFunction>,
    
    pub prim_num: usize, // contraction depth

    pub ele_offset: usize, // displacement from the start of Geometry.eles

    // easy indexing to AO representation
    pub num_ao: usize,
    pub ao_offset: usize
}

pub struct ContractedFunction {    
    // Angular momentum associate with this function
    pub lx: usize,
    pub ly: usize,
    pub lz: usize,

    // How to get the correct prim numbers from BasisSet
    pub coeff_offset: usize,
    pub exp_offset: usize,
}

impl BasisSet{
    pub fn normalize(&mut self) {
        //self.print_primtive_self_overlap("Primitive Self Overlap Before Normalization");
        //self.print_contracted_self_overlap("Contracted Self Overlap Before Normalization");
        self.get_prim_norms();
        self.normalize_primitives();

        //self.print_primtive_self_overlap("Primitive Self Overlap After Primitive Normalization");
        //self.print_contracted_self_overlap("Contracted Self Overlap After Primitive Normalization");
        self.get_contract_norms();
        self.normalize_contracted();

        //self.print_primtive_self_overlap("Primitive Self Overlap After Contracted Normalization");
        //self.print_contracted_self_overlap("Contracted Self Overlap After Contracted Normalization");
    }

    pub fn normalize_primitives(&mut self) {
        //println!("Primitive norms are: {:?}", prim_norms);
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {
                for i in 0..shell.prim_num {
                    let coeff_idx = primitive.coeff_offset + i;
                    let cur_coef = self.prim_coeffs[coeff_idx];
                    let cur_norm = self.prim_norms[coeff_idx];
                    self.prim_coeffs[coeff_idx] = cur_coef * cur_norm;
                }
            }
        }

        // Print Statements to Check Self-Overlap After Normalization
    }

    pub fn normalize_contracted(&mut self) {
        //println!("Contracted norms are: {:?}", contract_norms);
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {
                for i in 0..shell.prim_num {
                    let coeff_idx = primitive.coeff_offset + i;
                    let cur_coef = self.prim_coeffs[coeff_idx];
                    let cur_norm = self.contract_norms[coeff_idx];
                    self.prim_coeffs[coeff_idx] = cur_coef * cur_norm;
                }
            }
        }
    }

    pub fn get_prim_norms(&mut self) {
        /* From "Fundamentals of Molecular Integrals Evaluation" by
                Justin T. Fermann and Edward F. Valeev */
        let mut norms: Vec<f64> = vec![0.0; self.prim_coeffs.len()];
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {
                for i in 0..shell.prim_num {
                    let cur_coef = self.prim_coeffs[primitive.coeff_offset + i];
                    let cur_exp = self.prim_exp[primitive.exp_offset + i];
                    let cur_overlap = one_center_two_gaussian_integral(
                        primitive.lx, 
                        primitive.ly, 
                        primitive.lz, 
                        cur_exp, 
                        1.0, 
                        cur_exp, 
                        1.0);
                                        
                    self.prim_norms[primitive.coeff_offset + i] = sqrt(&(1.0 / cur_overlap));
                }
            }
        }
    }

    pub fn get_contract_norms(&mut self) {
        /* From "Fundamentals of Molecular Integrals Evaluation" by
                Justin T. Fermann and Edward F. Valeev */

        let mut norms: Vec<f64> = vec![0.0; self.prim_coeffs.len()];
        //println!("The Contracted Norm Contracted_Int_Sum are:");
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {

                let mut contracted_int_sum = 0.0;
                for i in 0..shell.prim_num {
                    let coeff1 = self.prim_coeffs[primitive.coeff_offset + i];
                    let exp1 = self.prim_exp[primitive.exp_offset + i];
                    
                    for j in 0..shell.prim_num {
                        let coeff2 = self.prim_coeffs[primitive.coeff_offset + j];
                        let exp2 = self.prim_exp[primitive.exp_offset + j];
                        contracted_int_sum += one_center_two_gaussian_integral(
                            primitive.lx,
                            primitive.ly, 
                            primitive.lz, 
                            exp1, 
                            coeff1, 
                            exp2, 
                            coeff2);
                    }
                }

                let prim_function_norm = 1.0 / sqrt(&(contracted_int_sum));

                for i in 0..shell.prim_num {
                    self.contract_norms[primitive.coeff_offset + i] = prim_function_norm;
                }
            }
        }
    }

    pub fn print(&self, mol: &molecule::Geometry) {
        /*Prints basis information for the current BasisSet, requires
        the Geometry to print the atom type */
        let header = "> Basis Sets <";
        println!("{:=^48}\n", header);
        for shell in self.shells.iter() {
            println!("{:?} Shell for atom: {}{}, l = {}",
            shell.gauss_type,
            molecule::ELEMENTS_TO_STR[&mol.eles[shell.ele_offset]],
            shell.ele_offset + 1,
            shell.l);

            for function in shell.functions.iter() {
                println!("{:2}lx: {}, ly: {}, lz: {}",
                "",
                function.lx, function.ly, function.lz);
                println!("{:4}Coefficients{:8}Exponents",
                "",
                "");
                for i in 0..shell.prim_num {
                    println!("{:6}{:>12.7} ({:^5}){:2}{:>12.7} ({:^5})",
                        "",
                        self.prim_coeffs[i + function.coeff_offset],
                        (i + function.coeff_offset).to_string(),
                        "",
                        self.prim_exp[i + function.exp_offset],
                        (i + function.exp_offset).to_string()
                    );
                }
            }
        }
        println!("{:=^48}\n", "");
    }

    fn print_primtive_self_overlap(&self, title: &str) {
        println!("\n{}", title);
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {
                for i in 0..shell.prim_num {
                    let coeff_idx = primitive.coeff_offset + i;
                    let exp_idx = primitive.exp_offset + i;

                    let cur_overlap = one_center_one_gaussian_integral(
                        primitive.lx, 
                        primitive.ly, 
                        primitive.lz, 
                        self.prim_exp[exp_idx], 
                        self.prim_coeffs[coeff_idx]);

                    println!("({:^5}) {:<12.8}", (coeff_idx + 1).to_string(), cur_overlap.to_string());
                }
            }
        }
    }

    pub fn print_contracted_self_overlap(&self, title: &str) {
        println!("\n{}", title);
        let mut count = 1;
        for shell in self.shells.iter() {
            for primitive in shell.functions.iter() {
                let mut contracted_int_sum = 0.0;
                for i in 0..shell.prim_num {
                    let coeff1 = self.prim_coeffs[primitive.coeff_offset + i];
                    let exp1 = self.prim_exp[primitive.exp_offset + i];
                    for j in 0..shell.prim_num {
                        let coeff2 = self.prim_coeffs[primitive.coeff_offset + j];
                        let exp2 = self.prim_exp[primitive.exp_offset + j];
                        contracted_int_sum += one_center_two_gaussian_integral(
                            primitive.lx,
                            primitive.ly, 
                            primitive.lz, 
                            exp1, 
                            coeff1, 
                            exp2, 
                            coeff2);
                    }
                }

                println!("({:^5}) {:<12.8}", count.to_string(), contracted_int_sum.to_string());
                count += 1;
            }
        }
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

pub fn load_basis (mol: &molecule::Geometry, basis_name: &str) -> BasisSet {
    /* This function is responsible for loading the basis information from the stored
    basis information as well as the molecular geometry to make sure that all the atoms
    have the apropriate cartesian gaussians
    
    TODO - adding ghost functions
    */

    let base_path = std::env!("CARGO_MANIFEST_DIR");
    let basis_path = &format!("{base_path}/basis/{basis_name}");
    let contents = &fs::read_to_string(basis_path).unwrap();

    let ele_basis_sets: parse_json::json_basis = serde_json::from_str(contents).unwrap();

    let mut mol_basis: BasisSet = build_mol_basis(&mol, &ele_basis_sets, BasisType::Cartesian);

    // Normalization - Does the primitive gaussians and then the contracted
    //mol_basis.normalize();

    // Testing to see if 1.18 pops up anywhere.. {
    mol_basis.get_prim_norms();
    mol_basis.normalize_primitives();

    println!("Before Contracted Normalization");
    mol_basis.print(&mol);
    mol_basis.get_contract_norms();
    mol_basis.normalize_contracted();

    println!("After Contracted Normalization");
    // } End Testing 
    mol_basis.print(&mol);

    // Recasting for immutable
    let mol_basis = mol_basis;

    return mol_basis;
}

pub fn build_mol_basis (geom: &molecule::Geometry, ele_basis_sets: &parse_json::json_basis, mode: BasisType) -> BasisSet {


    let func_type: &str = if mode == BasisType::Spherical {"gto_spherical"} else {"gto"};
    let mut coef_offset = 0 ;
    let mut exp_offset = 0 ;
    let mut ao_offset = 0 ;
    let mut mol_basis = BasisSet { 
        shells: Vec::new(),
        prim_coeffs: Vec::new(),
        prim_exp: Vec::new(),
        prim_norms: Vec::new(),
        contract_norms: Vec::new(),
        total_aos: 0 };
    
    for i in 0..geom.natoms {
        let cur_ele: u8 = geom.eles[i];
        let cur_shells: &Vec<parse_json::electron_shell> = &ele_basis_sets.elements[&(cur_ele as usize)].electron_shells;
        
        for shell in cur_shells.iter() {
            if shell.function_type == func_type {
                mol_basis.prim_exp.extend(shell.exponents.iter().copied());
                for ang in shell.angular_momentum.iter() {
                    let cur_basis_type: BasisType = if func_type == "gto_spherical" {BasisType::Spherical} else {BasisType::Cartesian};
                    
                    // before account for permutations of ang, this existed here
                    // mol_basis.prim_coeffs.extend(shell.coefficients[*ang].iter().copied());
                    
                    let mut contracted_functions: Vec<ContractedFunction> = Vec::new();

                    // Add inner loop here to account for each permutation of angular momentums for the basis sets
                    for lx in 0..(ang + 1) {
                        for ly in 0..(ang + 1 - lx) {
                            let lz = ang - lx - ly;
                            mol_basis.prim_coeffs.extend(shell.coefficients[*ang].iter().copied());
                            let cur_primitive: ContractedFunction = 
                                ContractedFunction { 
                                    lx: lx,
                                    ly: ly,
                                    lz: lz,

                                    coeff_offset: coef_offset,
                                    exp_offset: exp_offset
                                };
                            contracted_functions.push(cur_primitive);
                            coef_offset += shell.coefficients[*ang].len();
                        }
                    }

                    let mut cur_contract_shell: ContractedShell = 
                        ContractedShell {
                            gauss_type: cur_basis_type,
                            l: *ang,
                            prim_num: shell.exponents.len(),
                            functions: contracted_functions,
                            ele_offset: i,
                            num_ao: 0,
                            ao_offset: ao_offset
                        };
                    // deleted the times prim num to just get the number of contracted AOs present 
                    //let num_ao = cur_contract_shell.get_num_ao() * cur_contract_shell.prim_num;
                    cur_contract_shell.num_ao = cur_contract_shell.get_num_ao();

                    // tracking where in the ao block the shell is
                    ao_offset += cur_contract_shell.num_ao;
                    // immutable fixing
                    let cur_contract_shell = cur_contract_shell;
                    mol_basis.shells.push(cur_contract_shell);
                }
                exp_offset += shell.exponents.len();
            }
        }
    }

    mol_basis.total_aos = ao_offset;
    mol_basis.prim_norms = vec![1.0; mol_basis.prim_coeffs.len()];
    mol_basis.contract_norms = vec![1.0; mol_basis.prim_coeffs.len()];

    let mol_basis = mol_basis;
    return mol_basis;
}