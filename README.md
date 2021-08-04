# adaptation

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

#### HMC vs L-BFGS exam in Introduction (old test, no longer present in the paper, but the saved result of check_LBFGS.R is reused in the experiements in Section 3):
1. Optimization: posteriordb_check/phaseI_check/check_LBFGS.R
2. Stan phase I: posteriordb_check/phaseI_check/main.R
3. Plots: posteriordb_check/pic.R


#### Experiements in Section 3:
##### 100 repeats for experiments in Section 3
1. Stan phase I: posteriordb_check/phaseI_check/main_100_PhI.R
2. Pathfinder: posteriordb_check/phaseI_adapt_check/main_pf.R
3. ADVI: posteriordb_check/phaseI_ADVI/main_ADVI_100.R
##### Wasserstein distance computation and Plots for all experiements
4. Wasserstein distance: posteriordb_check/wasserstain_check.R
5. Plots: posteriordb_check/pic.R
##### Examples in Section 3.4
6. 100-dimensional Neal's Funnel: example/funnel100/funnel100.R
7. Ovarian example: example/ovarian/ovarian_test.R
8. Mcycle_gp_accel_gp: example/mcycle_gp/mcycle_gp.R
##### Birthday example in Section 3.5
9. example/Birthdays/model6_test.R


Authors
---------
| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Lu Zhang | lz2786@columbia.edu        | Ph.D.  Columbia Univeristy, Statistics |
| Bob Carpenter | bcarpenter@flatironinstitute.org       | Ph.D. Flatiron Institute |



Licensing
---------
* Code &copy; 2021, Lu Zhang, licensed under [BSD (3-clause)](https://opensource.org/licenses/BSD-3-Clause).

Notes
---------
* You are welcome to contact me for bug reporting.




