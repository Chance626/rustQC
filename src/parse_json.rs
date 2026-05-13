/*
Chance Brandt, The Ohio State University 2026

These are functions to help deserialize json files to load data into the runtime
*/

use serde::{Deserialize, Deserializer};
use std::collections::HashMap;

#[derive(Deserialize, Debug)]
pub struct json_basis {
    pub molssi_bse_schema: molssi_bse_schema,
    pub revision_description: String,
    pub revision_date: String,
    pub elements: HashMap<usize, elements>,
}

#[derive(Deserialize, Debug)]
struct molssi_bse_schema {
    pub schema_type: String,
    pub schema_version: String
}

#[derive(Deserialize, Debug)]
struct elements {
    pub electron_shells: Vec<electron_shell>,
    pub references: Vec<reference>,
}

#[derive(Deserialize, Debug)]
struct electron_shell {
    pub function_type: String,
    pub region: String,
    pub angular_momentum: Vec<usize>,

    #[serde(deserialize_with = "vec_str_to_f64")]
    pub exponents: Vec<f64>,

    #[serde(deserialize_with = "nested_vec_str_to_f64")]
    pub coefficients: Vec<Vec<f64>>,
}

#[derive(Deserialize, Debug)]
struct reference {
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

        for i in 1..self.elements.len() {
            println!("{}", i.to_string());
            println!("    electron_shells:");
            for shell in &self.elements[&i].electron_shells {
                
            }
        }
    }
}

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
