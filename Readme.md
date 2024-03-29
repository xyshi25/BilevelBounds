# Mixed Integer Bilevel Optimization with $k$-optimal Follower: A Hierarchy of Bounds

Xueyu Shi, Oleg A. Prokopyev, Ted K. Ralphs


## DESCRIPTION
The folder contains the code and dataset used in the paper "Mixed Integer Bilevel Optimization with $k$-optimal Follower: A Hierarchy of Bounds"


The paper in Optimization Online: http://www.optimization-online.org/DB_HTML/2020/06/7874.html


## SUPPORTED PLATFORMS
The code should be run in Linux platform

## BUILDING AND RUNING

### STEP 1
Extract the dataset through the command
```
tar -xzvf Data.tar.gz
```

### STEP 2
Modify the installation directionary of CPLEX in file "Makefile"

### STEP 3
Build the environment through the command
```
make all
```

### STEP 4
Run the experiments for the proposed methods through the command
```
chmod 777 run.sh
./run.sh
```

### STEP 5
Run the experiments for Mibs through the command (note that the solver Mibs should be installed first, see the instruction for Mibs in https://github.com/coin-or/MibS)
```
chmod 777 run_mibs.sh
./run_mibs.sh
```

### STEP 6
Clean the environment through the command
```
make clean
```


### STEP 7
Analyze the final computational results with the python scripts in foler "scripts". For example, to obtain the computational result in Tables 2 and 3, run
```
python3 scripts/anay_KIP.py
```
