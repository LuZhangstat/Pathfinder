# adaptation

This repository contains the code files of the analyses in "Pathfinder: Quasi-Newton posterior approximations alongoptimization paths"

Roadmap
---------
|Folder Name |     Intro            |
|:------ |:----------- |
|posteriordb| Subfolder posterior_database contains the database of package posteriordb|
|posteriordb_check| Codes for experiments before Section 3.6|
|example| Codes and data for Birthday example in Section 3.6 |
|utils| Utility functions |


Instructions
---------
#### HMC vs L-BFGS exam in Introduction:
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
##### Birthday example in Section 3.6
6. example/Birthdays/model6_test.R


Authors
---------
| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Lu Zhang | lz2786@columbia.edu        | Ph.D.  Columbia Univeristy, Statistics |
| Bob Carpenter | bcarpenter@flatironinstitute.org       | Ph.D., Flatiron Institute |



Licensing
---------
* Code &copy; 2021, Lu Zhang, licensed under [BSD (3-clause)](https://opensource.org/licenses/BSD-3-Clause).

Notes
---------
* You are welcome to contact me for bug reporting.




