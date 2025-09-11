#!/usr/bin/env python3
"""
Test the enhanced A-philic detector integration with app.py components
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_app_imports():
    """Test that app.py can import successfully with our enhanced detector"""
    print("Testing app.py imports...")
    
    try:
        # Test core imports from app.py
        from utils import parse_fasta, gc_content, reverse_complement, wrap
        from motifs import get_basic_stats
        from hyperscan_integration import all_motifs_refactored
        
        print("✅ Core app imports successful")
        return True
        
    except Exception as e:
        print(f"❌ App import failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_hyperscan_integration_aphilic():
    """Test that the enhanced A-philic detector works through hyperscan_integration"""
    print("\nTesting A-philic detection through hyperscan integration...")
    
    try:
        from hyperscan_integration import all_motifs_refactored
        
        # Test sequence with A-philic patterns
        test_seq = (
            "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
            "ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT"
            "GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA"
        )
        
        motifs = all_motifs_refactored(
            sequence=test_seq,
            sequence_name="test_sequence",
            nonoverlap=False,
            report_hotspots=True,
            calculate_conservation=False
        )
        
        print(f"Found {len(motifs)} total motifs")
        
        # Check for A-philic motifs
        aphilic_motifs = [m for m in motifs if m.get('Class') == 'A-philic DNA']
        print(f"Found {len(aphilic_motifs)} A-philic motifs")
        
        if aphilic_motifs:
            print("A-philic motif details:")
            for i, motif in enumerate(aphilic_motifs[:3]):  # Show first 3
                print(f"  Motif {i+1}:")
                for key in ['Start', 'End', 'Length', 'Normalized_Score', 'Subclass']:
                    if key in motif:
                        print(f"    {key}: {motif[key]}")
        
        return len(aphilic_motifs) > 0
        
    except Exception as e:
        print(f"❌ Hyperscan integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_motif_compatibility():
    """Test that A-philic motifs have all required fields for app.py"""
    print("\nTesting motif compatibility with app.py...")
    
    try:
        from hyperscan_integration import all_motifs_refactored
        
        test_seq = "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCATCGATCGCGCGCGCGATCG"
        
        motifs = all_motifs_refactored(
            sequence=test_seq,
            sequence_name="compatibility_test",
            nonoverlap=False,
            report_hotspots=True,
            calculate_conservation=False
        )
        
        # Check that all motifs have required fields
        required_fields = ['Sequence Name', 'Class', 'Subclass', 'Start', 'End', 'Length']
        
        for i, motif in enumerate(motifs):
            missing_fields = [field for field in required_fields if field not in motif]
            if missing_fields:
                print(f"❌ Motif {i} missing fields: {missing_fields}")
                return False
        
        print(f"✅ All {len(motifs)} motifs have required fields")
        
        # Test specific A-philic motifs
        aphilic_motifs = [m for m in motifs if m.get('Class') == 'A-philic DNA']
        
        if aphilic_motifs:
            # Check A-philic specific fields
            aphilic_fields = ['Normalized_Score', 'Actual_Score', 'Scoring_Method']
            sample_motif = aphilic_motifs[0]
            
            print("A-philic motif fields:")
            for field in aphilic_fields:
                value = sample_motif.get(field, 'MISSING')
                print(f"  {field}: {value}")
                
            # Check that the subclass is correctly set
            subclass = sample_motif.get('Subclass', '')
            if 'A-philic' in subclass:
                print("✅ A-philic subclass correctly set")
            else:
                print(f"⚠️ Unexpected subclass: {subclass}")
        
        return True
        
    except Exception as e:
        print(f"❌ Compatibility test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_basic_stats_integration():
    """Test that get_basic_stats works with A-philic motifs"""
    print("\nTesting basic stats integration...")
    
    try:
        from motifs import get_basic_stats
        from hyperscan_integration import all_motifs_refactored
        
        test_seq = "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCATCGATCGCGCGCGCGATCG"
        
        motifs = all_motifs_refactored(
            sequence=test_seq,
            sequence_name="stats_test",
            nonoverlap=False,
            report_hotspots=True,
            calculate_conservation=False
        )
        
        # Test basic stats calculation
        stats = get_basic_stats(test_seq, motifs)
        
        print("Basic stats:")
        for key, value in stats.items():
            print(f"  {key}: {value}")
        
        # Check that stats include expected fields
        expected_stats = ['Length', 'GC%', 'AT%', 'A', 'T', 'G', 'C', 'Motif Coverage %']
        missing_stats = [stat for stat in expected_stats if stat not in stats]
        
        if missing_stats:
            print(f"❌ Missing stats: {missing_stats}")
            return False
        
        print("✅ All expected stats present")
        return True
        
    except Exception as e:
        print(f"❌ Basic stats test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Enhanced A-philic App Integration Tests")
    print("=" * 50)
    
    tests = [
        test_app_imports,
        test_hyperscan_integration_aphilic,
        test_motif_compatibility,
        test_basic_stats_integration
    ]
    
    results = []
    for test in tests:
        try:
            results.append(test())
        except Exception as e:
            print(f"❌ Test {test.__name__} crashed: {e}")
            results.append(False)
    
    passed = sum(results)
    total = len(results)
    
    print(f"\n{'='*50}")
    print(f"Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("✅ All app integration tests passed!")
        print("The enhanced A-philic detector is ready for use in app.py")
    else:
        print("❌ Some integration tests failed!")
        for i, (test, result) in enumerate(zip(tests, results)):
            status = "✅" if result else "❌"
            print(f"  {status} {test.__name__}")
        sys.exit(1)