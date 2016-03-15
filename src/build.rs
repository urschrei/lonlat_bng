extern crate phf_codegen;

use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

const GENERATED_FILE: &'static str = "src/hexkeys.rs";

fn main() {
    let mut outfile = BufWriter::new(File::create(GENERATED_FILE).unwrap());
    write!(outfile, "static HEXKEYS: ::phf::Map<u32, &'static str> = ").unwrap();
    let mut builder = phf_codegen::Map::new();
    let rng = 000..1500;
    for val in rng {
        builder.entry(val, &format!("\"{:03x}\"", val));
    }
    builder.build(&mut outfile).unwrap();
    writeln!(outfile, ";").unwrap();
}
