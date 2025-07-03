use std::collections::HashMap;
use pyo3::{prelude::*, pyfunction};


// amino acid alphabeit
pub const AMINO_ACID_ALPHABET: &str = "ACDEFGHIKLMNPQRSTVWYXZ*";

// DNA + RNA alphabet
pub const XNA_ALPHABET: &str = "ACGTU";
// pub const DNA_ALPHABET: &str = "ACGT";
// pub const RNA_ALPHABET: &str = "ACGU";
pub const AMBIGUOUS_ALPHABET: &str = "ATGCUNRYSWKMBDHV";

// is amino acid string
#[pyfunction]
pub fn is_amino_acid_string(sequence: &str) -> bool {
    sequence.chars().all(|c| AMINO_ACID_ALPHABET.contains(c))
}

// is xna string
#[pyfunction]
pub fn is_xna_string(sequence: &str) -> bool {
    sequence.chars().all(|c| XNA_ALPHABET.contains(c))
}

// is ambiguous string
#[pyfunction]
pub fn is_ambiguous_string(sequence: &str) -> bool {
    sequence.chars().all(|c| AMBIGUOUS_ALPHABET.contains(c))
}


/// Complement a DNA base - using ambigous bases too
#[pyfunction]
pub fn complement_base(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        b'R' | b'r' => b'Y', // R = A or G; Y = C or T
        b'Y' | b'y' => b'R',
        b'M' | b'm' => b'K', // M = A or C; K = G or T
        b'K' | b'k' => b'M',
        b'S' | b's' => b'S', // S = G or C
        b'W' | b'w' => b'W', // W = A or T
        b'B' | b'b' => b'V', // B = C, G, or T; V = A, C, or G
        b'V' | b'v' => b'B',
        b'D' | b'd' => b'H', // D = A, G, or T; H = A, C, or T
        b'H' | b'h' => b'D',
        b'N' | b'n' => b'N', //  N  is wildcard, reverse of it is also wildcard
        _ => b'N' // Default to N for any other character
    }
}

/// Reverse complement a DNA sequence
#[pyfunction]
pub fn reverse_complement(sequence: &str) -> String {
    sequence.bytes()
        .rev()
        .map(complement_base)
        .map(|b| b as char)
        .collect()
}

/// Get reverse complement of a DNA sequence
#[pyfunction]
pub fn reverse_complement_seq(sequence: String) -> PyResult<String> {
    let seq = sequence.to_uppercase();
    // check all valid option (up to ambig)
    if !is_ambiguous_string(&seq) {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Sequence contains invalid nucleotides. Only A, T, C, G, U, R, Y, S, W, K, M, B, D, H, V, N are allowed."
        ));
    }
    
    Ok(reverse_complement(&seq))
}

/// Count nucleotides in a sequence
#[pyfunction]
pub fn count_nucleotides(sequence: String) -> PyResult<HashMap<char, usize>> {
    let mut counts = HashMap::new();
    
    for c in sequence.to_uppercase().chars() {
        if matches!(c, 'A' | 'T' | 'C' | 'G' | 'U' | 'N') {
            *counts.entry(c).or_insert(0) += 1;
        }
    }
    
    Ok(counts)
}

/// Calculate GC content of a sequence
#[pyfunction]
pub fn gc_content(sequence: String) -> PyResult<f64> {
    let seq = sequence.to_uppercase();
    let total = seq.len() as f64;
    let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
    
    if total == 0.0 {
        Ok(0.0)
    } else {
        Ok(gc_count / total)
    }
} 



