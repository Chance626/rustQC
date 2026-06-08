/*
Chance Brandt, The Ohio State 2026

Structure to store molecular information
*/

use std::collections::HashMap;
use std::sync::LazyLock;

use faer::traits::math_utils::sqrt;

pub static ELEMENTS: LazyLock<HashMap<&'static str, u8>> = LazyLock::new(|| {
    let mut m = HashMap::new();

    m.insert("H", 1);
    m.insert("He", 2);
    m.insert("Li", 3);
    m.insert("Be", 4);
    m.insert("B", 5);
    m.insert("C", 6);
    m.insert("N", 7);
    m.insert("O", 8);
    m.insert("F", 9);
    m.insert("Ne", 10);

    m.insert("Na", 11);
    m.insert("Mg", 12);
    m.insert("Al", 13);
    m.insert("Si", 14);
    m.insert("P", 15);
    m.insert("S", 16);
    m.insert("Cl", 17);
    m.insert("Ar", 18);
    m.insert("K", 19);
    m.insert("Ca", 20);

    m.insert("Sc", 21);
    m.insert("Ti", 22);
    m.insert("V", 23);
    m.insert("Cr", 24);
    m.insert("Mn", 25);
    m.insert("Fe", 26);
    m.insert("Co", 27);
    m.insert("Ni", 28);
    m.insert("Cu", 29);
    m.insert("Zn", 30);

    m.insert("Ga", 31);
    m.insert("Ge", 32);
    m.insert("As", 33);
    m.insert("Se", 34);
    m.insert("Br", 35);
    m.insert("Kr", 36);
    m.insert("Rb", 37);
    m.insert("Sr", 38);
    m.insert("Y", 39);
    m.insert("Zr", 40);

    m.insert("Nb", 41);
    m.insert("Mo", 42);
    m.insert("Tc", 43);
    m.insert("Ru", 44);
    m.insert("Rh", 45);
    m.insert("Pd", 46);
    m.insert("Ag", 47);
    m.insert("Cd", 48);
    m.insert("In", 49);
    m.insert("Sn", 50);

    m.insert("Sb", 51);
    m.insert("Te", 52);
    m.insert("I", 53);
    m.insert("Xe", 54);
    m.insert("Cs", 55);
    m.insert("Ba", 56);
    m.insert("La", 57);
    m.insert("Ce", 58);
    m.insert("Pr", 59);
    m.insert("Nd", 60);

    m
});

pub static ELEMENTS_TO_STR: LazyLock<HashMap<u8, &'static str>> = LazyLock::new(|| {
    ELEMENTS.iter().map(|(k, v)| (*v, *k)).collect()
});

pub struct Geometry {
    pub eles: Vec<u8>,
    pub coords: Vec<[f64; 3]>,
    pub natoms: usize,
    pub chrg: i8,
    pub mult: u8,
    pub nelec: u32
}

impl Geometry {
    pub fn print(&self) {
        let header = "> Molecular Geometry <";
        println!("{:=^48}\n", header);
        println!("charge: {}    multiplicity: {}\n", self.chrg, self.mult);
        for i in 0..self.natoms {
            println!("{:<5} {:3} {:12.8} {:12.8} {:12.8}", 
                    i + 1, 
                    ELEMENTS_TO_STR.get(&self.eles[i]).unwrap(), 
                    self.coords[i][0],
                    self.coords[i][1],
                    self.coords[i][2]);
        }
        println!("{:=^48}\n", "");
    }
}

/*
Geometry Util Functions
*/

#[inline]
pub fn get_cart_distance(loca: &[f64; 3], locb: &[f64; 3]) -> f64 {
    let mut square_sum = 0.0;
    for i in 0..3 {
        square_sum += (loca[i] - locb[i]).powi(2);
    }
    return sqrt(&square_sum);
}