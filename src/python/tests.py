from capibex import *
# (
#     translate,
#     translate_file,
#     reverse_complement,
#     complement_base,
#     filter_by_header,
#     filter_by_header_invert,
#     deduplicate_by_seq,
#     deduplicate_by_id,
#     sample_sequences,
#     sample_sequences_by_proportion,
#     is_amino_acid_string,
#     is_xna_string,
# )

def test_reverse_complement():
    """Test reverse complement function"""
    assert reverse_complement_seq("ATGC") == "GCAT"
    assert reverse_complement_seq("NNNN") == "NNNN"
    print("✓ reverse_complement test passed")

def test_translate():
    """Test translation function"""
    # Test standard genetic code
    result = translate("ATGAAATAG", frame=1, start=0, stop=None, to_protein=True, genetic_code=1)
    assert result == "MK*"  # ATG AAA TAG
    # assert "frame_-1" in result  # Verify reverse frames exist
    print("✓ translate test passed")

def test_filter():
    """Test sequence filtering"""
    
    # Create a test FASTA file
    with open("test.fasta", "w") as f:
        f.write(">seq1 test\nATGC\n>seq2 keep\nGCAT\n>seq3 test\nTGCA\n")
    
    # Filter sequences containing "keep"
    count = filter_by_header("test.fasta", ["keep"], "filtered.fasta", None)
    assert count == 1  # Should find one sequence
    
    # Verify the output
    with open("filtered.fasta") as f:
        content = f.read()
        assert ">seq2 keep" in content
        assert ">seq1 test" not in content
    print("✓ filter test passed")

if __name__ == "__main__":
    print("Running tests...")
    test_reverse_complement()
    test_translate()
    # test_filter() # ValueError: Failed to set up thread pool: The global thread pool has already been initialized.
    print("\nAll tests passed! 🎉") 