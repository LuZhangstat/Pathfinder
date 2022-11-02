# Pathfinder

This repository contains the code files of the analyses in "Pathfinder: Parallel quasi-Newton variational inference"

Roadmap
---------
|Folder Name |     Intro            |
|:------ |:----------- |
|posteriordb| Subfolder posterior_database contains the database of package posteriordb|
|posteriordb_check| Codes for experiments for examples from posteriordb |
|example| Codes and data for the Birthday example and other examples |
|utils| Utility functions |


Instructions
---------
(note: since the function **pdb_local()** requires folder posteriordb to be in a git repository, please use **git clone** to download the repository or configure a remote in the unzipped repository)

#### HMC vs L-BFGS exam in Introduction (old test, no longer present in the paper, but the result of check_LBFGS.R is reused in the experiments in Section 3):
1. Optimization: posteriordb_check/phaseI_check/check_LBFGS.R
2. Stan phase I: posteriordb_check/phaseI_check/main.R
3. Plots: posteriordb_check/pic.R


#### Experiements in Section 3:
##### 100 repeats for experiments in Section 3
1. Stan phase I: posteriordb_check/phaseI_check/main_100_PhI.R
2. Pathfinder: posteriordb_check/phaseI_adapt_check/main_pf.R
3. ADVI: posteriordb_check/phaseI_ADVI/main_ADVI_100.R
4. Laplace (in Appendix D): posteriordb_check/Laplace_check/main_laplace.R
##### Wasserstein distance computation and Plots for all experiements
5. Wasserstein distance: posteriordb_check/wasserstain_check.R
6. Plots: posteriordb_check/pic.R
##### Examples in Section 3.4
7. 100-dimensional Neal's Funnel: example/funnel100/funnel100.R
8. Ovarian example: example/ovarian/ovarian_test.R
9. Mcycle_gp_accel_gp: example/mcycle_gp/mcycle_gp.R
##### Birthday example in Section 3.5
10. example/Birthdays/model6_test.R


Authors
---------
| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Lu Zhang | lzhang63@usc.edu        | Ph.D.  University of Southern California, PPHS, Biostatistics |
| Bob Carpenter | bcarpenter@flatironinstitute.org       | Ph.D. Flatiron Institute |



Licensing
---------
* Code &copy; 2021, Lu Zhang, licensed under [BSD (3-clause)](https://opensource.org/licenses/BSD-3-Clause).

Notes
---------
* You are welcome to contact me for bug reporting.




