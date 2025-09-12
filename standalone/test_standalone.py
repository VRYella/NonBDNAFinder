#!/usr/bin/env python3
"""
Test script for the standalone NonBDNA Finder
"""

import os
import sys
import tempfile

# Add current directory to path
sys.path.insert(0, '.')

try:
    from nbdfinder_core import StandaloneNonBDNAFinder, analyze_fasta_file
    import pandas as pd
    print("✅ Successfully imported nbdfinder_core")
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)

def test_basic_functionality():
    """Test basic motif detection functionality"""
    print("\n🧪 Testing basic functionality...")
    
    # Create test sequence with known motifs
    test_fasta = """>test_g4_sequence
GGGTTTTGGGTTTTGGGTTTTGGGAAATTTCCCAAATTTCCCAAATTTCCCAAATTTCCC
>test_z_dna_sequence
CGCGCGCGCGCGCGCGCGCGCGCGCG
>test_a_tract_sequence
AAAAAAAAAAAAAAATGCGTAAAAAAAAAAAAA"""
    
    try:
        # Run analysis
        df, zip_data = analyze_fasta_file(test_fasta)
        
        print(f"   📊 Detected {len(df)} motifs")
        print(f"   📦 Zip package size: {len(zip_data):,} bytes")
        
        if len(df) > 0:
            class_counts = df['Class'].value_counts()
            print(f"   🏷️ Classes found: {dict(class_counts)}")
            
            # Check that we have the expected columns
            expected_cols = ['S.No', 'Sequence_Name', 'Class', 'Start', 'End', 'Length', 'Normalized_Score']
            missing_cols = [col for col in expected_cols if col not in df.columns]
            if missing_cols:
                print(f"   ⚠️ Missing columns: {missing_cols}")
            else:
                print(f"   ✅ All expected columns present")
            
            return True
        else:
            print(f"   ⚠️ No motifs detected")
            return False
            
    except Exception as e:
        print(f"   ❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_zip_contents():
    """Test that zip file contains expected files"""
    print("\n📦 Testing zip file contents...")
    
    test_fasta = """>simple_test
GGGTTTTGGGTTTTGGGTTTTGGG"""
    
    try:
        df, zip_data = analyze_fasta_file(test_fasta)
        
        # Extract zip to temporary directory
        import zipfile
        import io
        
        with tempfile.TemporaryDirectory() as temp_dir:
            zip_buffer = io.BytesIO(zip_data)
            with zipfile.ZipFile(zip_buffer, 'r') as zipf:
                zipf.extractall(temp_dir)
                
                # Check for expected files
                files = os.listdir(temp_dir)
                expected_files = ['nbdfinder_results.csv', 'nbdfinder_results.xlsx', 
                                'nbdfinder_results.gff3', 'nbdfinder_results_summary.txt']
                
                print(f"   📁 Files in zip: {files}")
                
                missing_files = [f for f in expected_files if f not in files]
                if missing_files:
                    print(f"   ⚠️ Missing files: {missing_files}")
                    return False
                else:
                    print(f"   ✅ All expected files present")
                    return True
                    
    except Exception as e:
        print(f"   ❌ Error testing zip contents: {e}")
        return False

def test_performance():
    """Test performance with a larger sequence"""
    print("\n⚡ Testing performance...")
    
    # Create a larger test sequence (1kb with repeating patterns)
    base_pattern = "GGGTTTTGGGTTTTGGGTTTTGGGAAACCCAAACCCAAACCCAAACCCATATATATATATATGCGCGCGCGCGCGC"
    large_sequence = base_pattern * 13  # ~1kb sequence
    
    test_fasta = f""">large_test_sequence
{large_sequence}"""
    
    try:
        import time
        start_time = time.time()
        
        df, zip_data = analyze_fasta_file(test_fasta)
        
        end_time = time.time()
        analysis_time = end_time - start_time
        
        print(f"   ⏱️ Analysis time: {analysis_time:.2f} seconds")
        print(f"   📏 Sequence length: {len(large_sequence):,} bp")
        print(f"   📊 Motifs detected: {len(df)}")
        print(f"   📦 Zip size: {len(zip_data):,} bytes")
        
        if analysis_time < 10:  # Should be much faster than 10 seconds
            print(f"   ✅ Performance acceptable")
            return True
        else:
            print(f"   ⚠️ Performance slower than expected")
            return False
            
    except Exception as e:
        print(f"   ❌ Performance test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("🧬 NonBDNA Finder Standalone - Test Suite")
    print("=" * 50)
    
    tests = [
        test_basic_functionality,
        test_zip_contents,
        test_performance
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print(f"\n📋 TEST SUMMARY")
    print(f"   Passed: {passed}/{total}")
    
    if passed == total:
        print(f"   ✅ All tests passed! Standalone version is working correctly.")
        return True
    else:
        print(f"   ❌ Some tests failed. Please check the issues above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)