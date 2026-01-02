# NonBDNAFinder Performance Benchmark Report
**Generated:** 2026-01-02 11:35:04
---

## Speed Performance

| Sequence Size | Time (s) | Speed (bp/s) | Motifs Detected | Density (motifs/kb) |
|---------------|----------|--------------|-----------------|---------------------|
| 1,000 bp | 0.20 | 4,889 | 54 | 54.00 |
| 5,000 bp | 0.42 | 11,910 | 209 | 41.80 |
| 10,000 bp | 0.72 | 13,972 | 412 | 41.20 |
| 50,000 bp | 3.89 | 12,870 | 2263 | 45.26 |
| 100,000 bp | 7.66 | 13,056 | 4499 | 44.99 |

**Average speed:** 11,340 bp/s

## Memory Usage

- Initial memory: 240.5 MB
- Peak memory: 253.3 MB
- Memory delta: 12.8 MB
- After cleanup: 243.6 MB
- Memory freed: 16.7 MB

## Accuracy Tests

| Test | Expected Class | Detected | Total Motifs | Status |
|------|----------------|----------|--------------|--------|
| Telomeric G4 | G-Quadruplex | 1 | 2 | ✓ PASS |
| HTT CAG repeat (normal) | Slipped-DNA | 0 | 1 | ✗ FAIL |
| HTT CAG repeat (expanded) | Slipped-DNA | 0 | 1 | ✗ FAIL |
| Z-DNA forming sequence | Z-DNA | 1 | 2 | ✓ PASS |
| A-tract curved DNA | Curved-DNA | 0 | 1 | ✗ FAIL |
| Palindrome (cruciform) | Cruciform | 1 | 2 | ✓ PASS |

**Success rate:** 50.0%

## Scalability

| Sequence Size | Time (s) | Speed (bp/s) | Memory (MB) | Status |
|---------------|----------|--------------|-------------|--------|
| 10,000 bp | 0.54 | 18,544 | 240.8 | ✓ SUCCESS |
| 50,000 bp | 2.97 | 16,847 | 240.2 | ✓ SUCCESS |
| 100,000 bp | 5.87 | 17,040 | 242.0 | ✓ SUCCESS |
| 500,000 bp | 29.09 | 17,188 | 250.9 | ✓ SUCCESS |

