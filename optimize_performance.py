#!/usr/bin/env python3
"""
Hyperscan Performance Optimizer for NBDFinder
=============================================

Pre-warms Hyperscan databases for faster startup and optimizes pattern compilation.
This should be run before starting the Streamlit app or API server for optimal performance.
"""

import sys
import os
import time

def prewarm_hyperscan_cache():
    """Pre-compile all Hyperscan databases for faster startup"""
    print("🔥 Pre-warming Hyperscan caches for optimal performance...")
    
    try:
        # Import all motif modules to trigger database compilation
        from motifs import (
            g_quadruplex, i_motif, z_dna, curved_dna, 
            slipped_dna, r_loop, cruciform_dna, triplex
        )
        
        # Test sequence to trigger all pattern compilations
        test_sequence = ("ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGG"
                        "CCCCCTCCCCCTCCCCCTCCCCATCGATCGCGCGCGCGATCGCACACACACAGCTGC"
                        "TGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTTGGGTTTAGGGGGGAGGGG")
        
        modules_to_prewarm = [
            ("G-Quadruplex", g_quadruplex),
            ("i-Motif", i_motif), 
            ("Z-DNA", z_dna),
            ("Curved DNA", curved_dna),
            ("Slipped DNA", slipped_dna),
            ("R-Loop", r_loop),
            ("Cruciform", cruciform_dna),
            ("Triplex", triplex)
        ]
        
        compiled_count = 0
        total_time = 0
        
        for name, module in modules_to_prewarm:
            start_time = time.time()
            
            try:
                # Run each module's detection to compile patterns
                if hasattr(module, 'find_g_quadruplex'):
                    module.find_g_quadruplex(test_sequence)
                elif hasattr(module, 'find_i_motif'):
                    module.find_i_motif(test_sequence)
                elif hasattr(module, 'find_z_dna'):
                    module.find_z_dna(test_sequence)
                elif hasattr(module, 'find_curved_dna'):
                    module.find_curved_dna(test_sequence)
                elif hasattr(module, 'find_slipped_dna'):
                    module.find_slipped_dna(test_sequence)
                elif hasattr(module, 'find_r_loop'):
                    module.find_r_loop(test_sequence)
                elif hasattr(module, 'find_cruciform_dna'):
                    module.find_cruciform_dna(test_sequence)
                elif hasattr(module, 'find_triplex'):
                    module.find_triplex(test_sequence)
                
                compile_time = time.time() - start_time
                total_time += compile_time
                compiled_count += 1
                
                print(f"   ✅ {name}: compiled in {compile_time:.3f}s")
                
            except Exception as e:
                print(f"   ⚠ {name}: compilation issue - {e}")
        
        # Get cache statistics
        try:
            from motifs.hyperscan_manager import get_hyperscan_cache_stats
            stats = get_hyperscan_cache_stats()
            
            print(f"\n📊 Pre-warming Results:")
            print(f"   • Modules processed: {compiled_count}/{len(modules_to_prewarm)}")
            print(f"   • Total compilation time: {total_time:.3f}s")
            print(f"   • Databases cached: {stats.get('database_cache_size', 0)}")
            print(f"   • Cache hits: {stats.get('cache_hits', 0)}")
            print(f"   ✅ Hyperscan optimization ready!")
            
        except ImportError:
            print(f"   ⚠ Cache statistics not available")
        
        return compiled_count > 0
        
    except Exception as e:
        print(f"   ❌ Pre-warming failed: {e}")
        return False

def optimize_for_production():
    """Apply production optimizations"""
    print("🚀 Applying production optimizations...")
    
    # Set environment variables for optimal performance
    os.environ['STREAMLIT_BROWSER_GATHER_USAGE_STATS'] = 'false'
    os.environ['STREAMLIT_BROWSER_COLLECTION_ENABLED'] = 'false'
    
    # Optimize Hyperscan settings
    prewarm_success = prewarm_hyperscan_cache()
    
    if prewarm_success:
        print("✅ NBDFinder optimized for production use!")
        print("🎯 Ready for fast Streamlit startup and API performance")
        return True
    else:
        print("⚠ Some optimizations failed - tool will still work but may be slower")
        return False

if __name__ == "__main__":
    print("🧬 NBDFinder Performance Optimizer")
    print("=" * 50)
    
    success = optimize_for_production()
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 Optimization complete! Start your application:")
        print("   • Streamlit: streamlit run app.py")
        print("   • REST API: python api.py")
    else:
        print("⚠ Optimization completed with warnings")
    
    sys.exit(0 if success else 1)