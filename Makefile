SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c

CP := g++

tfm_index_construct: seqana/lib/tfm_index_construct.cpp
	g++ $^ -l sdsl -l divsufsort -l divsufsort64 -o $@

