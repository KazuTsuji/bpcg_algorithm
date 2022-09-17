 # Files
 
 We include the following materials in the supplementary file. 

・FrankWolfe.jl 
You can reproduce the numerical experiments regarding convex optimization problems in the paper by              
importing this folder as a Julia package and implenting the "experiments.ipynb" in the "examples" folder. 
                                        
・ICML_kernel_herding_experiments
You can reproduce the experiments on the kernel herding methods.  The source codes are MATLAB files.

"bpcg_compair_methods.m" is the main file for comparing kernel quadrature methods. 
"time_bpcg_compair.m" is the file for comparing convergence of MMD for computational time. 

There are following correspondence between the file name and methods:

bcg_pairwise_linesearch.m・・・BPCG
bcg_pairwise_linesearch_lazified.m・・・BPCG_lazified
linesearch.m・・・kernel herding with linesearch
eqweight_herding.m・・・kernel herding with equal weights
PWF.m ・・・Pairwise CG for kernel herding
Away.m・・・Away CG for kernel herding
SBQ.m・・・Sequential Bayesian Quadrature
monte_carlo.m・・・Montecarlo
　








