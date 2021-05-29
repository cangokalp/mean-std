This repo consists of implementation of the proposed methods in "Mean-Standard Deviation Model For Minimum Cost Flow Problem" paper.

------

Networks used in this project can be found in the following repository - we use the benchmark suite NETGEN;

http://lime.cs.elte.hu/~kpeter/data/mcf/

------

usage: python solve_mcf.py [-e E] [-typ TYP [TYP ...]] [-exp EXP [EXP ...]]
                    [-t T [T ...]] [-l [L [L ...]]]

optional arguments:
  -h, --help          show this help message and exit
  -e E                experiment name
  -typ TYP [TYP ...]  list of network types to run the experiment on
  -exp EXP [EXP ...]  list of # of nodes, 2^k, input list of ks
  -t T [T ...]        list of instance extensions
  -l [L [L ...]]      weight parameter of the MSMCF instance

Example:
python solve_mcf.py -e graph_families -typ _sr_ -exp 12 13 -t a.min b.min 

This will run graph_family experiments on _sr_ network families for node sizes 2^12 and 2^13 on problem instances a and b

NOTE:

Note that the subproblems are solved using cplex, you would need to install cplex first. One can also change the solver option to another free solver in the code and solve using that.
