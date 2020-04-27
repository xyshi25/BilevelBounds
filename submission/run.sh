#!/bin/bash

# run experiments

# run experiments for the knapsack interdiction problem
mkdir -p Log/Log_KIP
./bin/KnapInterdiction

# run experiments for the bilevel vertex cover problem
mkdir -p Log/Log_BVC_Symm
mkdir -p Log/Log_BVC_ASymm
./bin/BVertexCover


# run experiments for the minimum spanning tree problem
mkdir -p Log/Log_BMST
./bin/BMST


# # to clean all the data
# make clean