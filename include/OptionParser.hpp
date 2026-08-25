#pragma once

#include "CommandOptions.hpp"

void help(char** argv, bool advanced=false);

// statistics about a GFA file
AppConfig main_stat(int argc, char** argv);
void help_stat(char** argv);

// Extract sequences along paths from a GFA
AppConfig main_seq(int argc, char** argv);
void help_seq(char** argv);

// extract sequences along paths from a GFA
AppConfig main_gfa2fa(int argc, char** argv);
void help_gfa2fa(char** argv);

// compute depth information from sequencing data for a GFA
AppConfig main_depth(int argc, char** argv);
void help_depth(char** argv);

// identify bubble in graph
AppConfig main_bubble(int argc, char** argv);
void help_bubble(char** argv, bool advanced=false);

// de-overlap graph
AppConfig main_deoverlap(int argc, char** argv);
void help_deoverlap(char** argv, bool advanced = false);

// collapse unitigs and bubbles (requires no-overlap graph, output by deoverlap)
AppConfig main_collapse(int argc, char** argv);
void help_collapse(char** argv, bool advanced=false);

// convert GFA/PAF to map format
AppConfig main_file2map(int argc, char** argv);
void help_file2map(char** argv);

// Liftover coordinates using a mapping file
AppConfig main_liftover(int argc, char** argv);
void help_liftover(char** argv, bool advanced=false);

// adjust mapping quality
AppConfig main_mapq_boost(int argc, char** argv);
void help_mapq_boost(char** argv, bool advanced=false);

// align reads to graph (requires no-overlap graph, output by collapse)
AppConfig main_align(int argc, char** argv);
void help_align(char** argv);

// restrict cuts
AppConfig main_res_cut(int argc, char** argv);
void help_res_cut(char** argv);

// split graph by connected components
AppConfig main_split(int argc, char** argv);
void help_split(char** argv);

// augment a GFA with VCF alleles
AppConfig main_augment(int argc, char** argv);
void help_augment(char** argv, bool advanced=false);

// remove weak UTG links using contig read components
AppConfig main_clean(int argc, char** argv);
void help_clean(char** argv);

// fill sample contig gaps from cross-sample paths
AppConfig main_gapfill(int argc, char** argv);
void help_gapfill(char** argv);

void finalize_opt_cfg(AppConfig& cfg);
