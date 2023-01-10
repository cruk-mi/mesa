# mesa 0.2.3

CHANGES
* `plotGeneHeatmap` now uses `bioMaRt` to find the gene details, based on the genome or mart supplied.
* `filterByOverlaps` now works when the qseaSet regions contains "chr" but the regions to filter by do not.
* Removed all references to `PooledControl` sample (previously suggested for normalisation).
* Various R CMD check warning issues fixed.
* `filterByOverlaps` no longer allows the option to take a list of window ID numbers.
* Fixed `addMedipsEnrichmentFactors` function inside `makeQset` when using the `"qseaPaired"` coverageMethod. 
* Fixed tests for `addMedipsEnrichmentFactors`

# mesa 0.2.2

BUG FIXES

* Update mixQsetSamples to match the mixArrayWithQset code (and now runs with latest version of R).
* Less stringent checking on the p values being zero bug in qsea, to minimise false positives.

KNOWN ISSUES
* `pull` is not working directly, `pullQset` must be used.

# mesa 0.2.1

Version set for initial release in v2.0.1 of the internal CBC Nextflow pipeline.
