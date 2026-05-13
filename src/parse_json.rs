/*
Chance Brandt, The Ohio State University 2026

These are functions to help deserialize json files to load data into the runtime
*/

use serde::{Deserialize, Deserializer};
use std::collections::HashMap;
use std::cmp::Ordering;

use crate::molecule;

#[derive(Deserialize, Debug)]
pub struct json_basis {
    pub molssi_bse_schema: molssi_bse_schema,
    pub revision_description: String,
    pub revision_date: String,
    pub elements: HashMap<usize, elements>,
}

#[derive(Deserialize, Debug)]
pub struct molssi_bse_schema {
    pub schema_type: String,
    pub schema_version: String
}

#[derive(Deserialize, Debug)]
pub struct elements {
    pub electron_shells: Vec<electron_shell>,
    pub references: Vec<reference>,
}

#[derive(Deserialize, Debug)]
pub struct electron_shell {
    pub function_type: String,
    pub region: String,
    pub angular_momentum: Vec<usize>,

    #[serde(deserialize_with = "vec_str_to_f64")]
    pub exponents: Vec<f64>,

    #[serde(deserialize_with = "nested_vec_str_to_f64")]
    pub coefficients: Vec<Vec<f64>>,
}

#[derive(Deserialize, Debug)]
pub struct reference {
    pub reference_description: String,
    pub reference_keys: Vec<String>
}

impl json_basis {
    pub fn print(&self) {
        println!("molssi_bse_schema:\n    schema_type: {}\n    schema_version: {}",
        self.molssi_bse_schema.schema_type,
        self.molssi_bse_schema.schema_version);
        println!("revision_description: {}",self.revision_description);
        println!("revision_date: {}", self.revision_date);
        println!("elements");

        for i in 1..(self.elements.len() + 1) {
            let num_to_ele = molecule::ELEMENTS_TO_STR.get(&(i as u8)).unwrap();
            println!("element: {}", num_to_ele.to_string());
            println!("    electron_shells:");

            for shell in &self.elements[&i].electron_shells {
                println!("        function_type: {}", shell.function_type);
                let ang_mom: String = shell.angular_momentum.iter()
                                    .map(|n| format!("{} ", n.to_string())).collect();
                println!("        angular_momentum: {}", ang_mom);
                let exp_str: String = shell.exponents.iter()
                                    .map(|n| format!("\n{:12}{}", "", n)).collect(); 
                println!("        exponents:{}", exp_str);
                if shell.angular_momentum.len()> 1 {
                    println!("        coefficients:");
                    for l in shell.angular_momentum.iter() {
                        let coef_str: String = shell.coefficients[0].iter()
                                        .map(|n| format!("\n{:12}{}", "", n)).collect();
                        println!("          for l = {}:{}", l, coef_str);
                    }
                } else {
                    let coef_str: String = shell.coefficients[0].iter()
                                    .map(|n| format!("\n{:12}{}", "", n)).collect();
                    println!("        coefficients:{}", coef_str);
                }
            }
            println!("    references:");
            for refer in self.elements[&i].references.iter() {
                println!("        {}", refer.reference_description);
                for refkey in refer.reference_keys.iter() {
                    println!("        - {}", refkey);
                }
            }
        }
    }
}


// Helper functions to parse the JSON formatted basis sets from BasisSetExchange

fn vec_str_to_f64<'de, D>(deserializer: D) -> Result<Vec<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let v = Vec::<String>::deserialize(deserializer)?;

    v.into_iter()
        .map(|s| s.parse::<f64>().map_err(serde::de::Error::custom))
        .collect()
}

fn nested_vec_str_to_f64<'de, D>(deserializer: D) -> Result<Vec<Vec<f64>>, D::Error>
where
    D: Deserializer<'de>
{
    let v = Vec::<Vec<String>>::deserialize(deserializer)?;

    v.into_iter().map(|inner| {
            inner.into_iter()
            .map(|s| s.parse::<f64>().map_err(serde::de::Error::custom)).collect()
        }).collect()
}

// ordering operators for the electron_shells
// made these for a print statement, but like... didn't need to...

impl PartialEq for electron_shell {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for electron_shell {}

impl PartialOrd for electron_shell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for electron_shell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.function_type.cmp(&other.function_type)
        .then(self.angular_momentum.cmp(&other.angular_momentum))
        .then_with(|| {
            self.exponents.len().cmp(&other.exponents.len())
            .then_with(|| {
                for i in 0..self.exponents.len() {
                    let exp_ord = self.exponents[i].total_cmp(&other.exponents[i]);
                    if exp_ord != Ordering::Equal {
                        return exp_ord;
                    }
                }
                Ordering::Equal
            })
        })
        .then_with(|| {
            self.coefficients.len().cmp(&other.coefficients.len())
            .then_with(|| {
                for i in 0..self.coefficients.len() {
                    let coef_len_ord = self.coefficients[i].len().cmp(&other.coefficients[i].len());
                    if coef_len_ord != Ordering::Equal {
                        return coef_len_ord;
                    }
                }
                Ordering::Equal
            })
            .then_with(|| {
                for i in 0..self.coefficients.len() {
                    for j in 0..self.coefficients[i].len() {
                        let coef_ord = self.coefficients[i][j].total_cmp(&other.coefficients[i][j]);
                        if coef_ord != Ordering::Equal {
                            return coef_ord;
                        }
                    }
                }
                Ordering::Equal
            })
        })
    }
}
