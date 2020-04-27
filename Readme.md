# A hierarchy of  bounds for bilevel mixed 0--1 programs


## DESCRIPTION
The folder contains the code and dataset used in the paper "A hierarchy of  bounds for bilevel mixed 0--1 programs"

The pdf version of the paper is named "Paper_bounds_for_bilevlMIP.pdf" 

## SUPPORTED PLATFORMS
The code should be run in Linux platform

## BUILDING AND INSTALLING

### STEP 1
Extract the dataset through command
```
	tar -xzvf Data.tar.gz
```

### STEP 2
Modify the installation directionary of CPLEX in file "Makefile"

### STEP 3
Build the environment through command
```
	make all
```

### STEP 4
Run the experiments for the proposed methods through the command
```
	chmod 775 run.sh
	./run.sh
```

### STEP 5
Run the experiments for Mibs through the command (note that the solver Mibs should be installed first, see the instruction for Mibs in https://github.com/coin-or/MibS)
```
	chmod 775 run_mibs.sh
	./run_mibs.sh
```

### STEP 6
Clean the environment through the command
```
	make clean
```