#!/usr/bin/env python3
"""
Test script for the domain scanning module.

This script validates the domain scanning functionality without requiring
actual InterProScan installation by using mock data.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd

# Add src to path
import sys
sys.path.insert(0, str(Path(__file__).parent / "src"))

from domain_scanning import DomainArchitectureScanner, DomainHit, DomainArchitecture


def create_mock_interproscan_output(output_file: Path) -> None:
    """Create mock InterProScan TSV output for testing"""
    mock_data = [
        "P12345\tMD5\t100\tPfam\tPF00001\tDomain1\t10\t50\t1e-10\tT\t2023-01-01\tIPR001234\tInterPro1\tGO:0001234\t",
        "P12345\tMD5\t100\tPfam\tPF00002\tDomain2\t60\t90\t1e-15\tT\t2023-01-01\tIPR005678\tInterPro2\tGO:0005678\t",
        "P67890\tMD5\t150\tPfam\tPF00001\tDomain1\t20\t60\t1e-12\tT\t2023-01-01\tIPR001234\tInterPro1\tGO:0001234\t",
        "P67890\tMD5\t150\tPfam\tPF00003\tDomain3\t80\t120\t1e-08\tT\t2023-01-01\tIPR009876\tInterPro3\tGO:0009876\t",
        "P11111\tMD5\t200\tPfam\tPF00001\tDomain1\t30\t70\t1e-20\tT\t2023-01-01\tIPR001234\tInterPro1\tGO:0001234\t"
    ]
    
    with open(output_file, 'w') as f:
        for line in mock_data:
            f.write(line + '\n')


def create_test_fasta(output_file: Path) -> None:
    """Create test FASTA file"""
    sequences = [
        ">P12345\nMETHIONINESTARTSEQUENCE",
        ">P67890\nANOTHERPROTEINSEQUENCE",
        ">P11111\nTHIRDPROTEINFORTESTING"
    ]
    
    with open(output_file, 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')


def test_domain_parsing():
    """Test InterProScan output parsing"""
    print("Testing domain parsing...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create mock InterProScan output
        tsv_file = temp_path / "test_interpro.tsv"
        create_mock_interproscan_output(tsv_file)
        
        # Create scanner with dummy InterProScan path
        scanner = DomainArchitectureScanner(
            interproscan_path=Path("/fake/interproscan"),
            num_cores=1
        )
        
        # Test parsing
        try:
            domain_map = scanner.parse_interpro_output(tsv_file)
            
            # Validate results
            assert 'P12345' in domain_map
            assert 'P67890' in domain_map
            assert 'P11111' in domain_map
            
            # Check P12345 domains
            p12345_domains = domain_map['P12345']
            assert 'PF00001' in p12345_domains
            assert 'PF00002' in p12345_domains
            assert 'PF00001,PF00002' in p12345_domains  # Supra-domain
            
            # Check P67890 domains
            p67890_domains = domain_map['P67890']
            assert 'PF00001' in p67890_domains
            assert 'PF00003' in p67890_domains
            assert 'PF00001,PF00003' in p67890_domains  # Supra-domain
            
            print("✓ Domain parsing test passed")
            print(f"  Found {len(domain_map)} proteins with domains")
            for protein, domains in domain_map.items():
                single_domains = [d for d in domains if ',' not in d]
                supra_domains = [d for d in domains if ',' in d]
                print(f"  {protein}: {len(single_domains)} single + {len(supra_domains)} supra-domains")
            
            return True
            
        except Exception as e:
            print(f"✗ Domain parsing test failed: {e}")
            return False


def test_supra_domain_generation():
    """Test supra-domain generation logic"""
    print("\nTesting supra-domain generation...")
    
    scanner = DomainArchitectureScanner(
        interproscan_path=Path("/fake/interproscan"),
        num_cores=1
    )
    
    # Test cases
    test_cases = [
        ([], []),  # Empty list
        (['PF00001'], []),  # Single domain
        (['PF00001', 'PF00002'], ['PF00001,PF00002']),  # Two domains
        (['PF00001', 'PF00002', 'PF00003'], ['PF00001,PF00002', 'PF00002,PF00003', 'PF00001,PF00002,PF00003']),  # Three domains
        (['PF00001', 'PF00002', 'PF00003', 'PF00004'], ['PF00001,PF00002', 'PF00002,PF00003', 'PF00003,PF00004', 'PF00001,PF00002,PF00003', 'PF00002,PF00003,PF00004'])  # Four domains (max 3 combo)
    ]
    
    all_passed = True
    for i, (input_domains, expected_supra) in enumerate(test_cases):
        try:
            result = scanner._generate_supra_domains(input_domains, max_length=3)
            if result == expected_supra:
                print(f"✓ Test case {i+1} passed: {input_domains} -> {result}")
            else:
                print(f"✗ Test case {i+1} failed: {input_domains}")
                print(f"  Expected: {expected_supra}")
                print(f"  Got: {result}")
                all_passed = False
        except Exception as e:
            print(f"✗ Test case {i+1} error: {e}")
            all_passed = False
    
    return all_passed


def test_sequence_preparation():
    """Test FASTA sequence preparation"""
    print("\nTesting sequence preparation...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create test FASTA
        input_fasta = temp_path / "test_input.fasta"
        create_test_fasta(input_fasta)
        
        output_fasta = temp_path / "prepared_sequences.fasta"
        
        scanner = DomainArchitectureScanner(
            interproscan_path=Path("/fake/interproscan"),
            num_cores=1
        )
        
        try:
            result_file = scanner.prepare_sequences_for_scanning([input_fasta], output_fasta)
            
            assert result_file.exists()
            assert result_file == output_fasta
            
            # Check content
            with open(result_file, 'r') as f:
                content = f.read()
                assert 'P12345' in content
                assert 'P67890' in content
                assert 'P11111' in content
            
            print("✓ Sequence preparation test passed")
            return True
            
        except Exception as e:
            print(f"✗ Sequence preparation test failed: {e}")
            return False


def test_sequence_validation():
    """Test sequence validation logic"""
    print("\nTesting sequence validation...")
    
    scanner = DomainArchitectureScanner(
        interproscan_path=Path("/fake/interproscan"),
        num_cores=1,
        min_sequence_length=10,
        max_sequence_length=100
    )
    
    # Mock sequence record for testing
    class MockRecord:
        def __init__(self, seq_str, seq_id="test"):
            self.seq = seq_str
            self.id = seq_id
    
    test_cases = [
        (MockRecord("ACDEFGHIKLMNPQRSTVWY"), True),  # Valid sequence
        (MockRecord("ACDEFG"), False),  # Too short
        (MockRecord("A" * 150), False),  # Too long
        (MockRecord("ACDEFGHIKLMNPQRSTVWYUOJ123"), False),  # Invalid characters
        (MockRecord("XXXXXXXXXXXXXXXXXXXXX"), False),  # Too many ambiguous
    ]
    
    all_passed = True
    for i, (record, expected) in enumerate(test_cases):
        try:
            result = scanner._validate_sequence(record, f"test_{i}")
            if result == expected:
                print(f"✓ Validation test {i+1} passed")
            else:
                print(f"✗ Validation test {i+1} failed: expected {expected}, got {result}")
                all_passed = False
        except Exception as e:
            print(f"✗ Validation test {i+1} error: {e}")
            all_passed = False
    
    return all_passed


def test_id_cleaning():
    """Test sequence ID cleaning"""
    print("\nTesting ID cleaning...")
    
    scanner = DomainArchitectureScanner(
        interproscan_path=Path("/fake/interproscan"),
        num_cores=1
    )
    
    test_cases = [
        ("sp|P12345|PROTEIN_HUMAN", "P12345"),
        ("tr|Q98765|PROTEIN_MOUSE", "Q98765"),
        ("simple_id", "simple_id"),
        ("id with spaces", "id_with_spaces"),
        ("123invalid", "seq_123invalid"),
        ("normal|pipe|format", "normal")
    ]
    
    all_passed = True
    for input_id, expected in test_cases:
        try:
            result = scanner._clean_sequence_id(input_id)
            if result == expected:
                print(f"✓ ID cleaning test passed: '{input_id}' -> '{result}'")
            else:
                print(f"✗ ID cleaning test failed: '{input_id}' -> '{result}' (expected '{expected}')")
                all_passed = False
        except Exception as e:
            print(f"✗ ID cleaning test error: {e}")
            all_passed = False
    
    return all_passed


def test_domain_architecture_creation():
    """Test DomainArchitecture object creation"""
    print("\nTesting domain architecture creation...")
    
    scanner = DomainArchitectureScanner(
        interproscan_path=Path("/fake/interproscan"),
        num_cores=1
    )
    
    # Create test domain hits
    hits = {
        'P12345': [
            DomainHit('P12345', 'PF00001', 10, 50, 1e-10),
            DomainHit('P12345', 'PF00002', 60, 90, 1e-15)
        ]
    }
    
    try:
        architectures = scanner.create_domain_architecture_objects(hits)
        
        assert 'P12345' in architectures
        arch = architectures['P12345']
        
        assert isinstance(arch, DomainArchitecture)
        assert arch.protein_id == 'P12345'
        assert arch.single_domains == ['PF00001', 'PF00002']
        assert arch.supra_domains == ['PF00001,PF00002']
        assert arch.domain_positions['PF00001'] == (10, 50)
        assert arch.domain_positions['PF00002'] == (60, 90)
        assert arch.total_length == 90
        
        print("✓ Domain architecture creation test passed")
        return True
        
    except Exception as e:
        print(f"✗ Domain architecture creation test failed: {e}")
        return False


def main():
    """Run all tests"""
    print("=== Domain Scanning Module Tests ===\n")
    
    tests = [
        test_supra_domain_generation,
        test_id_cleaning,
        test_sequence_validation,
        test_sequence_preparation,
        test_domain_parsing,
        test_domain_architecture_creation
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"✗ Test {test_func.__name__} crashed: {e}")
            results.append(False)
    
    # Summary
    print("\n=== Test Summary ===")
    passed = sum(results)
    total = len(results)
    
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("✓ All tests passed!")
        return 0
    else:
        print(f"✗ {total - passed} tests failed")
        return 1


if __name__ == "__main__":
    exit(main())