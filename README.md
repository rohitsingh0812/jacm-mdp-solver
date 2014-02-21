jacm-mdp-solver
===============

To run experiments:

$ bash experiments.sh <low> <high>

e.g.

$ bash experiments.sh 3 10

This will create a lot of folders in this directory with results. Every folder will have a strategy.

For each n value between <low> and <high>, to experiments are done:

1. with equal probability p=0.05 for each process to get hungry (results in n_1 folder)
2. with increasing probabilities from 0.05 -> 0.12 -> 0.19 -> ... (results in n_2 folder)

output of the program is stored in ./n_m/output.txt (m=1 or 2)
The last row of the output has state-size, value and time taken (and some other details)
