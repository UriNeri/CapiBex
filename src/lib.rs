use pyo3::prelude::*;
use pyo3::Python;

// Import our modules
mod translate;
mod dedup;
mod filter;
mod sample;
mod seq_utils;
mod fastx;

// Expose the PyO3 modules
#[pymodule]
#[pyo3(name = "capibex")]
fn capibex(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(translate::translate, m)?)?;
    m.add_function(wrap_pyfunction!(translate::translate_file, m)?)?;


    m.add_function(wrap_pyfunction!(dedup::deduplicate_by_seq, m)?)?;
    m.add_function(wrap_pyfunction!(dedup::deduplicate_by_id, m)?)?;
    
    m.add_function(wrap_pyfunction!(filter::filter_by_header, m)?)?;
    m.add_function(wrap_pyfunction!(filter::filter_by_header_invert, m)?)?;

    m.add_function(wrap_pyfunction!(sample::sample_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(sample::sample_sequences_by_proportion, m)?)?;

    // Add direct functions
    
    m.add_function(wrap_pyfunction!(seq_utils::reverse_complement_seq, m)?)?;
    m.add_function(wrap_pyfunction!(seq_utils::count_nucleotides, m)?)?;
    m.add_function(wrap_pyfunction!(seq_utils::gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(translate::list_genetic_codes, m)?)?;

    m.add_function(wrap_pyfunction!(seq_utils::is_amino_acid_string, m)?)?;
    m.add_function(wrap_pyfunction!(seq_utils::is_xna_string, m)?)?;
    m.add_function(wrap_pyfunction!(seq_utils::is_ambiguous_string, m)?)?;
    m.add_function(wrap_pyfunction!(seq_utils::complement_base, m)?)?;

    m.add_class::<fastx::PyFastxReader>()?;
    m.add_class::<fastx::Record>()?;
    m.add_wrapped(wrap_pyfunction!(fastx::py_parse_fastx_file))?;
    m.add_wrapped(wrap_pyfunction!(fastx::parse_fastx_string))?;
    m.add_wrapped(wrap_pyfunction!(fastx::normalize_seq))?;
    m.add_wrapped(wrap_pyfunction!(fastx::reverse_complement))?;
    m.add_wrapped(wrap_pyfunction!(fastx::py_decode_phred))?;
    m.add("FastxParseError", py.get_type::<FastxParseError>())?;

    Ok(())
}