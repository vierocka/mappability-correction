# Tests: roadmap

| status | test | note |
|---|---|---|
| done | trimming impact (subsets) | small; paper keeps untrimmed |
| done | personal MANE impact (subsets) | fcM required for HLA |
| done | filtering impact vs MS (full data) | no AUC gain from ×funique |
| done | inflation tail (full data) | correction enriched in extremes, mild |
| done | HLA typing + reconcile (subsets) | 12/12 |
| done | end-to-end: calls → personal MANE (3 lines) | byte-identical to published |
| to do | full server run on one library (HPC) | user runs; commands in ../ROADMAP.md |
| to do | trimming impact on full libraries | subsets only so far |
| to do | pseudogene-decoy simulation | design in ../pseudogene_problem/decoys.md |
| to do | spliced-read funique vs exon-flank funique | junction reads missing in the exon-flank design |
| to do | stranded counting (`-s 2`) | antisense overlaps in the APM locus |
| to do | new read length (e.g. L88 for trimmed JY) | needs the pseudogene-evidence table at that length |
