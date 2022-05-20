# Tunneling and sequence analysis
This directory contains code to demonstrate the use of tunneling within the field of sequence analysis.

## Requirements
- GNU make
- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [gurobi-optimizer](https://www.gurobi.com/products/gurobi-optimizer/)

## What is contained
- A bundle of programs for the construction of a tunneled FM-index using ILP solver Gurobi
- A bundle of programs handling the recovery of the original string from the tunneled FM-index

## Program compilation
To compile all programs, just call `make`.

The software uses the Succinct Data Structure Library [sdsl-lite](https://github.com/simongog/sdsl-lite) by Simon Gog.

After installing sdsl lite, you can either
- set up an environment variable called `SDSLLITE` that specifies the path to the library
  using the command `export SDSLLITE="path/to/sdsllite/"`.
- modify the file `Make.helper` in the root directory such that the uppermost path points to the Make.helper file
  of your sdsl-lite installation.

## ILP-reduction and tunneled FM-index construction

### Installation
The program `tfm_index_construct.x` uses the ILP solver Gurobi.
Install it in the directory BWT-Tunneling-ILP, following the instruction at https://www.gurobi.com/documentation/quickstart.html.

After installing the ILP solver and obtaining the license, set up environment variables (GUROBI_HOME, PATH, LD_LIBRARY_PATH, GRB_LICENSE_FILE) using the commands in the mentioned documentation

### What is contained
- A program `tfm_index_construct.x` to construct a tunneled FM-index from a file.
  The file must not contain nullbytes, further information can be found by executing the program without parameter.
- A program `tfm_index_invert.x` which can be used to recover the original string from which the FM-index was built.
  Further information can be found by executing the program without parameter.
