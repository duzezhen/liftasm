# Changelog

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
