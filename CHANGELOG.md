# Changelog

## v0.1.3-r31 (2026-08-29)

### Added

- Added `--min_freq` and `--max_occ` to filter frequent minimizers and speed up alignment.
- Added `--min_ali_len` to filter short, potentially false alignments.

## v0.1.3-r30 (2026-08-28)

### Added

- Added VCF augmentation for renamed CTG inputs, overlapping variants, and sample-specific variant paths.
- Added `MT`/`TF` bubble annotations and primary haplotype mosaic-focused graphs.

### Fixed

- Fixed stale repeat-path edges and boundary-node inflation of homologous alignment conflict groups.

### Improved

- Improved parallel scheduling to avoid FIFO head-of-line blocking.

## v0.1.3-r29 (2026-08-27)

### Improved

- Updated automatic node-support thresholds to 2 for 3-4 haplotypes and 3 otherwise.
- Split fully low-support contigs that share nodes with a confirmed misassembly.

## v0.1.3-r28 (2026-08-26)

### Changed

- Increased the default `collapse --path_sim` from `0.5` to `0.9`.

## v0.1.3-r27 (2026-08-25)

### Fixed

- Fixed `--min_eq` filtering after alignments were split at graph node boundaries.
- Fixed VCF haplotype assignment to follow GFA P-lines when available.

## v0.1.3-r26 (2026-08-24)

### Improved

- Simplified terminal misassembly detection with per-component haplotype support and cross-component confirmation.
- Added final P-line-guided repeat-bubble normalization for contig collapse.

## v0.1.3-r25 (2026-08-23)

### Changed

- `collapse` now reads UTG/CTG inputs and optional VCFs from a header-driven configuration table.

### Improved

- Improved P-line/SN-aware bubble and homologous-path enumeration.
- Added clear errors for undefined segments referenced by GFA P/L/A lines.

## v0.1.3-r24 (2026-08-20)

### Fixed

- Fixed incorrect contig layout during gap filling.

## v0.1.3-r23 (2026-08-19)

### Added

- Added source-aware node-level gap walks with sample-haplotype soft masking.
- Added `--full_phase` for phase-supported haplotype exchange at overlapping diploid gaps.
- Added separate haplotype-1 and haplotype-2 bubble VCFs on the final primary contig coordinates.

### Changed

- Use spanning contigs only to prove connections and choose shared-node boundaries; fill sequence now comes from the graph walk.

## v0.1.3-r22 (2026-08-16)

### Improved

- Improved gap-fill boundary alignment, overlap splicing, and coordinate reporting.

## v0.1.3-r21 (2026-08-12)

### Fixed

- Improved overlapping-contig cleanup when some contigs covered by filled gaps were not removed.

## v0.1.3-r20 (2026-08-11)

### Improved

- Improved phase-aware gap selection and exact boundary refinement.
- Relocated terminal misassemblies before final gap filling.
- Fixed oversized gap bridges and incomplete multi-haplotype splitting.

## v0.1.3-r19 (2026-08-07)

### Improved

- Improved gap-fill boundary filtering, phasing, and debug output.
- Simplified gap-fill options and overlap filtering.
- Improved contig-mode collapse output and naming.

## v0.1.3-r18 (2026-08-06)

### Improved

- Improved phase-aware gap selection and topology-guided local boundary alignment.
- Added exact base-pair gap splicing with repeat-aware fallback to shared graph nodes.
- Improved phased primary output and removal of redundant small components and contigs.

## v0.1.3-r17 (2026-08-04)

### Fixed

- Prevented misassemblies and repeats from joining unrelated components.
- Fixed strand-insensitive clustering of bubble paths.
- Improved phase-aware gap filling to reduce switch errors.

## v0.1.3-r16 (2026-08-03)

### Added

- Added the `gapfill` subcommand to fill contig gaps using phased paths from other samples.

## v0.1.3-r15 (2026-07-28)

### Added

- Added contig mode to directly collapse haplotype 1/2 `p_ctg` GFA files from multiple samples.

## v0.1.3-r14 (2026-07-19)

### Added

- Added `--trim_len` and `--trim_ovlp` to trim small alignment overlaps.

### Changed

- Used `hi >= params_.min_similarity` to evaluate homologous path support in the first step to include the dead end path.

## v0.1.3-r13 (2026-07-17)

### Improved

- Reduced homologous path enumeration complexity with `open_walk()`.

## v0.1.3-r12 (2026-07-17)

### Fixed

- Fixed several cycle-inducing issues in `collapse`.

## v0.1.3-r11 (2026-06-14)

### Added

- Fallback to `align_short_mm2_()` when `align_mm2_()` returns no hits.

### Improved

- Reduced the time complexity of `build_bubble_branch_pairs_()`.

## v0.1.3-r10 (2026-06-12)

### Added

- Added `--min_ali_ratio` and `--min_mapq` parameters to the `deoverlap` and `collapse` subcommands to filter low-confidence alignments that may introduce cycles.
- Added the `liftasm split` subcommand to separate GFA file(s) by connected component.

### Changed

- Prevented merging of uniquely connected nodes when their `SN:Z:` sample origins differ.
- Simplified homologous path/bubble comparison by using the longest path as the reference and aligning all other paths against it.

## v0.1.3-r9 (2026-06-04)

### Added

- Added `-n, --name` to the `bubble`, `deoverlap`, and `collapse` subcommands.
- Merged segments now record sample origins using the `SN:Z:` tag.
- `--min_eq` in `collapse` now accepts multiple values for iterative collapse.
-  Added `-x, --preset` to the `deoverlap` and `collapse` subcommands to specify the alignment mode.

### Improved

- Faster path traversal using sample-origin information.
- VCF output from `bubble` now includes per-sample genotypes.

## v0.1.3-r8 (2026-05-29)

### Changed

- Removed the `--min_frac` parameter from the `mapq_boost` subcommand.

## v0.1.3-r7 (2026-05-26)

### Improved

- Added multithreading support for cut point building, propagation, and pruning (`build_propagate_prune_cuts_`).
- Reduced the time complexity of cut point propagation.

## v0.1.3-r6 (2025-05-25)

### Added

- Added `--abnormal_cut` parameter to filter consecutive short cut points after propagation, reducing cycle formation in the graph.

## v0.1.3-r5 (2026-05-24)

### Changed

- Updated `enumerate_paths_greedy_DFS_()` to prevent the same node from being used more than once within a single DFS path, regardless of orientation. Previously, `node+` and `node-` were treated as different vertices; now they are treated as the same underlying node during path-cycle detection.

## v0.1.3-r4 (2026-05-23)

### Added

- Added `--repeat_mask` parameter to the `collapse` subcommand to mask repetitive sequences within nodes, reducing the generation of cycles in the graph.

### Notes

- Enabling `--repeat_mask` increases alignment runtime because repetitive regions must be masked and corresponding CIGAR strings modified during each CIGAR generation step (~1.54× slower).

## v0.1.3-r3 (2026-05-20)

### Improved

- Reduced the time complexity of `has_conflicting_pair_rule_()`.
- Updated `build_rulemap_()` to use multithreading based on node connectivity.

### Fixed

- Fixed an issue where standalone nodes with no incoming or outgoing edges were incorrectly removed during `expand_and_rewire_edges_()` and `remove_unused_nodes_()`.

## v0.1.3-r2 (2025-05-18)

### Fixed

- Added a counter to `enumerate_paths_greedy_DFS_()` to prevent greedy DFS from getting trapped in local paths.
- Added `has_conflicting_pair_rule_()` to detect conflicting pair rules and prevent loops caused by shifted alignment windows in repeated node-pair alignments.
- Updated `build_rules_from_pair_windows_()` to support multiple mappings for the same key.
