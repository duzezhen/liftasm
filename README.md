# liftasm

## Introduction

Alignment-based coordinate liftover for assemblies and GFA graphs.

## Requirements

Please note the following requirements before building and running the software:

- `Linux` operating system
- GCC with `C++20` support (GCC ≥ 8)
- `CMake ≥ 3.20`
- `zlib`, `bzip2`, `liblzma` (xz)
- `Git` (for submodules)

Optional:
- `libcurl` (`htslib` will be built without `libcurl` support if not available)

## Installation

<!-- ### Install via conda

```shell
conda create -n syndiv
conda activate syndiv
# Install SynDiv with all dependencies
conda install -c bioconda -c conda-forge -c duzezhen syndiv
``` -->

### Building on Linux

Use the following script to build the software:

```shell
git clone https://github.com/duzezhen/liftasm.git
cd liftasm
git submodule update --init --recursive
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 5
./build/liftasm
```

If you would like to enable CPU-specific optimizations (e.g., `-march=native`), you can compile with:

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLIFTASM_ENABLE_NATIVE=ON
```

## Usage

### Input Files

<!-- To quickly get started, you will need two input files: `aligns` files and `syri.out` files. Once you have obtained these files, make sure to prepare the Reference genome and configuration file.

* Reference Genome
* configuration file

Please note that the chromosome names in the query genome must match those in the reference genome.

Sample names and file paths should not contain special characters, especially `_`.

```shell
# configuration file
sample1 sample1.fa sample1.aligns sample1.syri.out
sample2 sample2.fa sample2.aligns sample2.syri.out
...
sampleN sampleN.fa sampleN.aligns sampleN.syri.out
```

[configuration_url]: https://github.com/JiaoLab2021/SynDiv/wiki/Configuration-file

File should be separated by tabs. The code examples for generating `aligns` and `syri.out` files can be found on the [wiki][configuration_url].

### Running

Before running the software, it is recommended to set the maximum number of open files using the `ulimit -n <number>` command. The maximum number of open files can be calculated based on the number of genomes (`n`) and the number of threads (`t`) using the following formula:

```shell
number = 2nt
```

For convenience, let's assume the following file names for the input:

* `refgenome.fa`
* `configuration.txt`

**One-Click Generation**

```shell
ulimit -n 50000
SynDiv -r refgenome.fa -c configuration.txt &
```

**Note**

[Manual-execution_url]: https://github.com/JiaoLab2021/SynDiv/wiki/Manual-execution

See the [wiki][Manual-execution_url] for step-by-step manual execution of SynDiv and calucate Syn-Fst. -->

## License

MIT