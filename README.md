# liftasm

`liftasm` is a toolkit for de-overlapping and collapsing [hifiasm](https://github.com/chhylp123/hifiasm) assembly graphs. It works directly with `hifiasm` GFA output and can merge multiple GFA graphs, especially assemblies from different tissues of the same individual.

The central `collapse` command has two modes:

- **p_utg mode:** accepts one or more `p_utg` GFAs, de-overlaps them when necessary, and collapses homologous graph sequence.
- **p_ctg mode:** accepts paired haplotype-contig (`p_ctg`) GFAs from different tissues or samples and merges them into a combined graph. The final graph can be passed to `gapfill` to fill supported gaps and identify variants between haplotypes or tissues.

liftasm also writes coordinate maps that can be used for `liftover` and `mapq_boost`.

## Installation

liftasm requires Linux, GCC 12.2 or newer, CMake 3.20 or newer, zlib, bzip2, liblzma, and POSIX threads. `libcurl` is optional.

```bash
git clone https://github.com/duzezhen/liftasm.git
cd liftasm
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
```

The executable is `build/liftasm`. On a dedicated machine, add `-DLIFTASM_ENABLE_NATIVE=ON` to the CMake command to enable CPU-specific optimization.
Use `-DLIFTASM_ENABLE_CURL=OFF` to disable HTSlib remote URL access explicitly.

### Bash completion

To enable Bash completion, add the following lines to `~/.bashrc`. Replace `/path/to/liftasm` with the absolute path to this repository:

```bash
export PATH="/path/to/liftasm/build:$PATH"

# liftasm Bash completion
_liftasm_completion_file="/path/to/liftasm/completion/liftasm.bash"
if [[ -f "$_liftasm_completion_file" ]]; then
    source "$_liftasm_completion_file"
fi
unset _liftasm_completion_file
```

Reload the configuration with `source ~/.bashrc`. You can then type `liftasm` followed by <kbd>Tab</kbd> to complete subcommands.

## Quick start: collapse a p_utg graph

Use **UTG mode** when the input is one or more `p_utg`-style GFA graphs. You can pass the original overlap graph directly: `collapse` de-overlaps it internally when necessary. Running `deoverlap` first is optional, not required.

```bash
liftasm collapse \
  -t 16 \
  -g blood.p_utg.gfa \
  -n blood
```

To merge graphs from different tissues of one individual, list GFAs and names in the same order:

```bash
liftasm collapse \
  -t 16 \
  -g blood.p_utg.gfa heart.p_utg.gfa \
  -n blood heart
```

The two names after `-n` correspond to the two GFAs in the same order.

## Quick start: collapse haplotype contigs, then fill gaps

Use **CTG mode** when each tissue or sample has paired haplotype-contig GFAs. `--c1` contains haplotype 1 GFAs, `--c2` contains haplotype 2 GFAs, and the entries must be in the same order. The `-n` list labels those pairs.

```bash
liftasm collapse \
  -t 16 \
  --c1 blood.hap1.p_ctg.gfa heart.hap1.p_ctg.gfa \
  --c2 blood.hap2.p_ctg.gfa heart.hap2.p_ctg.gfa \
  -n blood heart
```

Run gap filling on the final CTG-mode graph:

```bash
liftasm gapfill \
  -t 16 \
  -g out.iter3.collapse.gfa
```

`gapfill` writes haplotype-specific GFAs for each input tissue or sample, primary haplotype GFAs, primary-coordinate bubble VCFs for hap1 and hap2, and reports named `out.gapfill.tsv`, `out.gapfill.relocations.tsv`, and `out.gapfill.html`.

## Main commands

### `collapse`

This is the central command. It searches graph bubbles and longer homologous paths, then merges supported equivalent sequence into shared nodes.

The most useful controls are `--iterations` (how many rounds to run), `--min_match` and `--min_ali_ratio` (how strict alignments must be), `--min_mapq` (minimum alignment confidence), `--min_eq` (minimum exact-match length used to create segment cut points), `--depth` and `--paths` (how broadly bubbles are explored), and `--ctg_min_len` (ignore short contigs in CTG mode). The default `--min_eq 1000,3,3` uses a conservative first round, then allows shorter exact matches in later rounds to refine the graph.

### `gapfill`

This command is intended for the final graph from CTG-mode `collapse`. A real spanning contig first establishes the connection and trusted shared-node boundaries; the gap sequence is then chosen by a source-aware node walk rather than copied from that contig. `--full_phase` lists samples with chromosome-scale phase information. Their haplotypes may be paired across a gap using local bubble phase, and only listed full-phase samples can support those connections. Other samples are joined within each haplotype.

### `deoverlap`

Use `deoverlap` when you need a non-overlap GFA as a standalone output, for another tool, or for inspection. It is not a prerequisite for `collapse -g`.

```bash
liftasm deoverlap \
  -t 16 \
  -g blood.p_utg.gfa \
  -n blood
```

### `liftover` and `file2map`

`liftover` moves BED intervals through one or more map files. Maps can come from `collapse`, `deoverlap`, or `file2map`.

```bash
liftasm liftover \
  -t 8 \
  -m out.iter1.collapse.map out.iter2.collapse.map out.iter3.collapse.map \
  -b regions.bed \
  -o regions.lifted.bed
```

To lift between assemblies or to a reference, first convert a `minimap2` PAF. Generate PAF with `-c` so that it contains CIGAR information.

```bash
minimap2 -cx asm5 --secondary=no -o aln.paf reference.fa assembly.fa
liftasm file2map -i aln.paf -o aln.paf.map

liftasm liftover \
  -t 8 \
  -m graph.map aln.paf.map \
  --paf aln.paf \
  -b regions.bed -o regions.reference.bed
```

### `mapq_boost`

`mapq_boost` is used when reads are aligned to their own diploid genome and homologous sequence between haplotypes causes MAPQ to fall to 0. It uses a homology map to correct this homology-driven ambiguity. The input must contain all secondary alignments for each read, either as separate secondary records or in alternative-alignment tags such as BWA's `XA`. Alignments for the same read must be adjacent; use `samtools sort -n` if necessary.

```bash
samtools sort -n -@ 8 -o reads.namesort.bam reads.bam
liftasm mapq_boost \
  -t 8 --io_threads 8 --name_check \
  -m out.iter1.collapse.map out.iter2.collapse.map out.iter3.collapse.map \
  -i reads.namesort.bam \
  -o reads.boosted.bam
```

It only changes alignments at or below `--mapq_low` (default: 3), never raises them beyond `--mapq_cap` (default: 60), and leaves unsupported ambiguous mappings unchanged.

## Further reading

[`docs/USAGE.md`](docs/USAGE.md) provides a detailed introduction to every `liftasm` subcommand and its parameters.

## Limitations

The complexity of a collapsed graph depends strongly on the original assembly graph. When multiple `p_utg` graphs are merged, some regions may form tangles. Improving these complex regions is ongoing work.

## License

MIT
