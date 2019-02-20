#[macro_use]
extern crate clap;
extern crate debruijn;
extern crate dna_io;
extern crate fnv;

use clap::{App};
use std::collections::{HashMap, HashSet};
use std::hash::BuildHasherDefault;
use std::io::BufRead;

use std::error::Error;

use std::cmp::min;

use std::fs::File;
use std::io::{BufReader, Read};

use debruijn::*;
use debruijn::kmer::*;
use debruijn::dna_string::*;
use fnv::FnvHasher;
type FnvHashMap<V,T> = HashMap<V,T, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<V> = HashSet<V, BuildHasherDefault<FnvHasher>>;

static mut KMER_SIZE: usize = 0;

fn main() {
    let parameters = load_params();
    let (het_kmers, pairs) = load_kmers(&parameters);
    identify_het_kmers(het_kmers, pairs, &parameters);
}


fn identify_het_kmers(het_kmers: FnvHashSet<u64>, pairs: Vec<(u64,u64)>, params: &Params) {
    let mut primary_het_kmers: FnvHashSet<u64> = FnvHashSet::default();
    let mut secondary_het_kmers: FnvHashSet<u64> = FnvHashSet::default();


    let reader = dna_io::DnaReader::from_path(&params.primary);
    for record in reader {
        for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
            let krc = k.rc();
            let to_hash = min( k.to_u64(), krc.to_u64() );
            if het_kmers.contains(&to_hash) {
                primary_het_kmers.insert(to_hash);
            }
        }
    }
    let reader = dna_io::DnaReader::from_path(&params.secondary);
    for record in reader {
        for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
            let krc = k.rc();
            let to_hash = min( k.to_u64(), krc.to_u64() );
            if het_kmers.contains(&to_hash) {
                secondary_het_kmers.insert(to_hash);
            }
        }
    }
    let mut stat_count = [[[[0u64; 2]; 2]; 2]; 2];
    for (het1, het2) in &pairs {
        let in_prim1: usize = match primary_het_kmers.contains(&het1) {true => 1, false => 0,};
        let in_prim2: usize = match primary_het_kmers.contains(&het2) {true => 1, false => 0,};
        let in_sec1: usize = match secondary_het_kmers.contains(&het1) {true => 1, false => 0,};
        let in_sec2: usize = match secondary_het_kmers.contains(&het2) {true => 1, false => 0,};
        stat_count[in_prim1][in_prim2][in_sec1][in_sec2] += 1;
    }
    println!("Error cases");
    println!("{} with both alts in primary and both alts in secondary",stat_count[1][1][1][1]);
    println!("{} with both alts in primary and one in secondary", 
        stat_count[1][1][1][0]+stat_count[1][1][0][1]);
    println!("{} with both alts in primary and neither in secondary", stat_count[1][1][0][0]);
    println!("{} with one alt in primary and the same alt in secondary, other alt not seen", 
        stat_count[1][0][1][0]+stat_count[0][1][0][1]);
    println!("{} with both alts in secondary and one alt in primary", 
        stat_count[1][0][1][1]+stat_count[0][1][1][1]);
    println!("{} with neither alt in primary and both alts in secondary",stat_count[0][0][1][1]);
    println!("{} with neither alt in primary and one alt in secondary",
        stat_count[0][0][0][1]+stat_count[0][0][1][0]);
    println!("{} with neither alt in primary or secondary",stat_count[0][0][0][0]);
    println!("total error cases = {}",
        stat_count[1][1][1][1] +
        stat_count[1][1][1][0] +
        stat_count[1][1][0][1] +
        stat_count[1][1][0][0] + 
        stat_count[1][0][1][0] +
        stat_count[0][1][0][1] +
        stat_count[1][0][1][1] +
        stat_count[0][1][1][1] +
        stat_count[0][0][1][1] +
        stat_count[0][0][0][1] +
        stat_count[0][0][1][0] + 
        stat_count[0][0][0][0]);
    println!();
    println!("True unphased case");
    println!("{} with one alt in primary and neither in secondary = true unphased",stat_count[0][1][0][0]+stat_count[1][0][0][0]);
    println!();
    println!("True phased case");
    println!("{} with one alt in primary and the other alt in secondary",stat_count[0][1][1][0]+stat_count[1][0][0][1]);

    
}

fn load_kmers(params: &Params) -> (FnvHashSet<u64>, Vec<(u64,u64)>){
    let mut kmers: FnvHashSet<u64> = FnvHashSet::default();
    let mut pairs: Vec<(u64,u64)> = Vec::new();
    let f = File::open(&params.het_kmers).expect("Unable to open het kmer file");
    let f = BufReader::new(f);
    for line in f.lines() {
        let line = line.expect("Unable to read line");
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let first = DnaString::from_dna_string(&tokens[0]);
        let second = DnaString::from_dna_string(&tokens[2]);
        unsafe {
            if KMER_SIZE == 0 && first.len() != 0 && first.len() == second.len() {
                KMER_SIZE = first.len();
            } else if first.len() != KMER_SIZE || second.len() != KMER_SIZE {
                panic!("kmer sizes must be consistent");
            }
        }
        let to_hash1: KmerX = first.get_kmer(0);
        let to_hash2: KmerX = second.get_kmer(0);
        kmers.insert(to_hash1.to_u64());
        kmers.insert(to_hash2.to_u64());
        pairs.push((to_hash1.to_u64(), to_hash2.to_u64()));
    }
    (kmers, pairs)
}



type KmerX = VarIntKmer<u64, KX>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct KX;

impl KmerSize for KX {
    fn K() -> usize {
        unsafe {
            KMER_SIZE
        }
    }
}

#[derive(Clone)]
struct Params {
    primary: String,
    secondary: String,
    het_kmers: String,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let fasta = params.value_of("primary_fasta").unwrap();
    let secondary = params.value_of("secondary_fasta").unwrap(); 
    let het_kmers = params.value_of("het_kmers").unwrap();
    Params{
        primary: fasta.to_string(),
        secondary: secondary.to_string(),
        het_kmers: het_kmers.to_string(),
    }
}
