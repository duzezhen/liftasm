# Using liftasm

This guide explains how to choose a command, what each command produces, and which parameters are worth changing in practice. It deliberately focuses on decisions a user needs to make; use `liftasm <command> --help` or `--advanced` for the complete interface and exact defaults.

## Inputs and threads

Graph commands take more than filenames because liftasm needs to know which sample each graph belongs to.

- `collapse -i FILE` reads one sample per row from a header-driven, whitespace-delimited table. Column order is arbitrary. The header selects UTG or CTG mode, and sample names must be unique.
- Other graph commands that accept multiple GFAs use `-n, --name` for input labels in the same order as their GFA files.
- `-t, --threads` controls compute threads. Use it with `collapse`, `deoverlap`, `gapfill`, `liftover`, and `mapq_boost`. `collapse` may need about 5 GB RAM for each thread, so choose a lower value on memory-limited machines.

## Which workflow fits the input?

| Input you have | Command to start with | What it does |
| --- | --- | --- |
| One or more hifiasm `p_utg` GFAs | `collapse -i samples.utg.tsv` | UTG mode: directly collapses homologous graph sequence. Overlap is handled internally when needed. |
| Paired hap1/hap2 `p_ctg` GFAs from different tissues or samples | `collapse -i samples.ctg.tsv` | CTG mode: compares paired contigs and builds a combined graph. |
| Final graph from CTG mode | `gapfill -g` | Uses spanning-contig evidence and source-aware graph walks to fill conservative contig gaps. |
| BED coordinates plus graph/assembly relationships | `liftover -m -b` | Transfers BED intervals through one or more coordinate maps. |
| SAM/BAM with low MAPQ in known homologous regions | `mapq_boost -m` | Removes homology-only ambiguity from MAPQ decisions. |
| A standalone non-overlap graph is needed | `deoverlap -g` | Converts an overlap GFA and writes a map. |

## `collapse`: the main graph operation

`collapse` finds equivalent sequence in bubbles and longer homologous paths, then replaces supported duplicate sequence with shared graph nodes. It creates a `.map` file at every collapse iteration; these maps are the link between the simplified graph and the original coordinates.

### Input table format

The first non-empty, non-comment line is a required header. Column names are lowercase and may appear in any order. Spaces and tabs are accepted as separators, although tabs are recommended. Blank lines and lines beginning with `#` are ignored. Because whitespace separates fields, sample names and paths cannot themselves contain spaces or tabs.

| Column | Required in | Description |
| --- | --- | --- |
| `sample` | UTG and CTG | Unique sample name. |
| `gfa` | UTG | Unitig GFA; cannot be combined with CTG columns. |
| `hap1_gfa` | CTG | Haplotype-1 contig GFA. |
| `hap2_gfa` | CTG | Haplotype-2 contig GFA; must occur together with `hap1_gfa`. |
| `vcf` | Optional CTG | Per-sample VCF. Use `.` when the sample has no VCF. |

Column order is arbitrary, but only these schemas are accepted:

- UTG: `sample`, `gfa`
- CTG: `sample`, `hap1_gfa`, `hap2_gfa`
- CTG with variants: `sample`, `hap1_gfa`, `hap2_gfa`, `vcf`

Unknown or repeated header columns, incompatible column combinations, duplicate sample names, and data rows whose field count differs from the header terminate with exit code 1. The error identifies the configuration file, line number, and reason; row errors also show the sample when available.

### UTG mode: input `p_utg` GFA directly

Use a table containing `sample` and `gfa` for one or more unitig graphs. A separate `deoverlap` command is not needed: `collapse` accepts an overlap or a non-overlap GFA and de-overlaps internally when required.

```text
sample    gfa
blood     blood.p_utg.gfa
heart     heart.p_utg.gfa
```

```bash
liftasm collapse \
  -t 16 \
  -i samples.utg.tsv
```

With the default output prefix and three rounds, the final files include:

- `out.iter3.collapse.gfa`: graph with sequence.
- `out.iter3.collapse.noseq.gfa`: graph topology without sequence.
- `out.iter3.collapse.map`: coordinate mapping for downstream use.

### CTG mode: paired haplotype-contig inputs

Use CTG mode for different tissues from one individual, or for multiple samples, with one haplotype-1 and one haplotype-2 `p_ctg` GFA on each input row.

```text
sample    hap1_gfa                     hap2_gfa
blood     blood.hap1.p_ctg.gfa         blood.hap2.p_ctg.gfa
heart     heart.hap1.p_ctg.gfa         heart.hap2.p_ctg.gfa
```

```bash
liftasm collapse \
  -t 16 \
  -i samples.ctg.tsv
```

This command first writes collapsed GFAs for each tissue, then writes the combined graph (`out.iter3.collapse.gfa` by default). CTG mode also saves raw and filtered component alignments, which is useful when checking how contigs were connected. To carry variants into the output, add the `vcf` header column and place each VCF on the same row as its sample; use `.` when a sample has no VCF.

### Parameters users commonly change

| Option | What it controls | When to change it |
| --- | --- | --- |
| `--iterations INT` `[3]` | Number of collapse rounds. | Increase only when later rounds still make useful merges. For comma-separated parameter lists, the last value is reused if the list is shorter than the number of rounds. |
| `--min_match FLOAT[,..]` `[0.9]` | Required matching-base fraction in a retained alignment. | Raise for stricter merges; lower only for more divergent assemblies. |
| `--min_ali_ratio FLOAT[,..]` `[0.05,0.01,0.001]` | Minimum aligned fraction of a sequence. | Raise to ignore short local similarities. |
| `--min_mapq INT` `[30]` | Minimum mapping quality for graph alignments. | Raise when repetitive regions introduce weak matches. |
| `--min_jaccard FLOAT[,..]` `[0.8]` | How similar two path memberships must be before grouping. | Raise to make path grouping more conservative. |
| `--min_eq INT[,..]` `[1000,3,3]` | Minimum exact-match length used to introduce graph cut points. | The default uses `1000` bp in round 1 to avoid over-cutting, then `3` bp in later rounds to refine existing homologous regions. Raise it for fewer, more conservative cuts; lower it only when finer collapse is needed and the graph is well supported. |
| `--depth INT` `[100kb]` | Maximum graph distance explored inside a bubble. | Increase for long bubbles; it can increase run time. |
| `--paths INT` `[0]` | Maximum alternative paths kept per bubble. | Raise only if a real bubble has more alternatives than the default retains. |
| `--ctg_min_len INT` `[2Mb]` | Minimum input contig length in CTG mode. | Lower when short contigs are important; `0` disables filtering. |
| `--anchor_only` | Align only shared-read anchors in CTG mode. | Use when full-contig alignment is unnecessary or too expensive. |

For difficult repetitive graphs, `--repeat_mask MIN_LEN,MAX_PERIOD,MAX_MISMATCH` controls tandem-repeat detection. UTG mode masks repeat matches during alignment; CTG mode uses the final contig P-lines to normalize repeat-only bubbles. `--no_cx` disables complex-region filtering and can permit unreliable complex-graph merges.

## `gapfill`: use source-aware node walks to close contig gaps

`gapfill` takes the final output of CTG-mode `collapse`, not the original `p_ctg` GFAs. A real spanning contig must support both sides of a gap and establish trusted shared-node boundaries. The sequence between those anchors is then selected by a source-aware graph walk: it prefers one sample/haplotype source and softly avoids nodes already used by the other haplotype, while shared homologous nodes remain available. Ambiguous connections are rejected.

```bash
liftasm gapfill \
  -t 16 \
  -g out.iter3.collapse.gfa
```

With the default prefix, the report is `out.gapfill.tsv`. The command also writes a relocation report and an HTML report, plus per-tissue, primary, and low-quality haplotype GFAs. Bubbles that map unambiguously to the final primary node paths are written separately as `out.primary.hap1.gapfill.bubbles.vcf` and `out.primary.hap2.gapfill.bubbles.vcf`; each file uses its own primary haplotype as the reference coordinate system and standard 1-based VCF coordinates.

`--full_phase` accepts one or more base sample names (the first column of the `collapse` input table). Full-phase targets may exchange haplotype contigs according to local bubble-phase evidence, and only listed full-phase samples may supply their bridges. Samples not listed are joined within the original haplotype.

| Option | What it controls | When to change it |
| --- | --- | --- |
| `--full_phase STR ...` | Samples assembled with chromosome-scale phase information. | List Hi-C, Pore-C, or trio-phased samples; omit ordinary locally phased assemblies. |
| `--depth INT` `[100k]` | Maximum graph distance searched for one bubble or gap walk. | Increase only when a valid walk exceeds the default node depth. |
| `--DFS_guard INT` `[1M]` | Maximum visited states for one bubble or gap walk. | Increase when a complex but valid region reaches the search limit. |
| `--min_match FLOAT` `[0.9]` | Minimum matching-base fraction in relocation and deduplication alignments. | Raise it for more conservative sequence checks. |
| `--min_ali_ratio FLOAT` `[0.05]` | Minimum aligned fraction in relocation and deduplication alignments. | Raise it to ignore short local hits. |
| `--ms_len INT` `[10M]` | Minimum low-support sequence at a contig end. | Lower it to inspect shorter terminal errors. |
| `--ms_sim FLOAT` `[0.7]` | Fraction that one other component must cover. | Raise it for more conservative splitting. |
| `--ms_haps INT\|auto` `[auto]` | Minimum supporting haplotypes, computed per component. | Use an integer to replace the automatic 2/3-haplotype rule. |
| `--min_mapq INT` `[30]` | Minimum mapping quality in relocation and deduplication alignments. | Raise it when repeats create ambiguous hits. |
| `--min_contig INT` `[1Mb]` | Shortest contig considered as evidence. | Raise it to avoid noisy short contigs. |
| `--max_gap INT` `[10Mb]` | Largest graph gap eligible for filling. | Lower it when only short joins are trusted. |
| `--max_overlap FLOAT` `[0.1]` | Largest allowed overlap between the two target contigs. | Lower it to reject heavily overlapping target pairs. |
| `--min_overlap INT` `[1Mb]` | Required shared-node overlap between a bridge and each target contig. | Raise for stricter connection evidence. |
| `--min_similarity FLOAT` `[0.7]` | Required shared-node similarity for homologous spanning evidence. | Raise for closer samples or stricter sibling-path evidence. |
| `--phase_len INT` `[10]` | Use only bubbles whose every internal path is at most this many bp. | Normally leave unchanged. |
| `--phase_win` | Local window used once for bubble phase evidence. | Defaults to 500 kb; a boundary with no bubble evidence is reported as `.`. |
| `--dedup_sim FLOAT` `[0.95]` | Similarity threshold for removing redundant contigs/components. | Lower only if legitimate near-duplicate contigs must be retained. |

## `deoverlap`: create a standalone non-overlap graph

`deoverlap` is useful when another program needs a non-overlap graph or when you want to inspect this intermediate representation. It is optional before UTG-mode `collapse`.

```bash
liftasm deoverlap \
  -t 16 \
  -g blood.p_utg.gfa \
  -n blood
```

It produces `.deoverlap.gfa`, `.deoverlap.noseq.gfa`, and `.deoverlap.map` files. `--min_match` `[0.9]`, `--min_ali_ratio` `[0.001]`, and `--min_mapq` `[30]` filter alignments used to resolve overlaps. Higher values make transformations more conservative. `--trim_len` and `--trim_ovlp` are advanced controls for trimming small overlaps between long alignments.

## `file2map` and `liftover`: move BED coordinates

A `.map` contains pairs of equivalent coordinate intervals. `collapse` and `deoverlap` write maps automatically. `file2map` creates one from GFA A-lines or a PAF alignment, allowing the graph map to be linked to contigs or a reference.

```bash
minimap2 -cx asm5 --secondary=no reference.fa assembly.fa > assembly.paf
liftasm file2map -i assembly.paf -o assembly.paf.map

liftasm liftover \
  -t 8 \
  -m out.iter3.collapse.map assembly.paf.map \
  --paf assembly.paf \
  -b intervals.bed \
  -o intervals.reference.bed
```

PAF conversion requires minimap2 `-c` so each alignment has a `cg` CIGAR tag. The reference and query names in that PAF must be different, and cannot contain `:`, `-`, or `+`.

| Option | What it controls |
| --- | --- |
| `liftover -m FILE ...` | One or more maps, in any order that creates a connected coordinate path. |
| `liftover --min_frac FLOAT` `[0.9]` | Minimum fraction of a BED interval that must map. |
| `liftover --cm_max_hops INT` `[8]` | Maximum number of map links followed. Lower it to restrict indirect mappings. |
| `liftover --regex STR` | Keep only destination labels matching a regular expression. |
| `liftover --check -r FASTA` | Check sequence consistency; useful for validating important mappings. |
| `file2map --min_len INT` `[50kb]` | Minimum PAF alignment length retained. |
| `file2map --min_mapq INT` `[5]` | Minimum PAF MAPQ retained. |

## `mapq_boost`: correct homology-driven MAPQ loss

`mapq_boost` is intended for reads aligned to their own diploid genome, where homologous sequence between haplotypes can reduce MAPQ to 0. It changes MAPQ only when the alternative placements of a read correspond through the supplied homology map; unsupported multi-mapping reads remain unchanged.

All secondary alignments for each read must be visible to the command, either as separate secondary records or in alternative-alignment tags such as BWA's `XA`. Records with the same read name must be consecutive. If the input BAM is coordinate-sorted, sort it by query name first.

```bash
samtools sort -n -@ 8 -o reads.namesort.bam reads.bam
liftasm mapq_boost \
  -t 8 --io_threads 8 --name_check \
  -i reads.namesort.bam \
  -o reads.boosted.bam \
  -m out.iter3.collapse.map
```

| Option | What it controls |
| --- | --- |
| `--mapq_low INT` `[3]` | Only reads with MAPQ at or below this value are candidates for boosting. |
| `--mapq_cap INT` `[60]` | Hard ceiling for a boosted MAPQ. |
| `--name_check` | Verify that all records for each read are adjacent. |
| `--sub_ovlp_frac FLOAT` `[0.95]` | How much two alignments must overlap on the read before they are compared together. |
| `--cm_max_hops INT` `[8]` | Search depth through the homology map. |
| `--io_threads INT` `[4]` | HTSlib input/output threads; independent of compute `-t`. |
| `-g FASTA -e MOTIF [MOTIF ...] --rc` | Optional restriction-site evidence for Hi-C, Pore-C, or CiFi data. `^` marks the cut position. |

The advanced `--K_*`, `--W_*`, and `--close_as_eps` options change how MAPQ, alignment score, match length, edit distance, and restriction-site evidence are weighted. Keep their defaults unless you have benchmark data for a different scoring model.

## Checking a run

For any command, start by saving the invoked command and inspecting the generated files:

```bash
liftasm collapse --help
liftasm collapse --advanced
liftasm gapfill --help
```

For graph operations, keep the final `.gfa`, `.map`, and gapfill reports together. For liftover, use `--check` on consequential intervals. For MAPQ boosting, inspect the change count in the program log and compare downstream results with the unmodified BAM.
