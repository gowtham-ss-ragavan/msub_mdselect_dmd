# Mean Subtraction and Mode Selection in DMD

### Getting started
- Save all .m in *examples* as .nb files within the same directory

#### Setup for analyzing cavity flow data:
- Data in a specific folder
  	- Create a folder named *data* in the same directory as *examples*
  	- Create two folders within *data*: *raws* and *pruned*
  	- Put all 4 files named [*Cavity_Rexxk.mat*](https://ucsb.app.box.com/s/7t44ootww9b3axe3oqyngttokl7f75lr) from [Hassan's cavity mixing paper](https://arxiv.org/abs/1903.10044) inside the folder *raws*
- Install MATLAB and MATLink. This will allow *Mathematica* to process the *.mat* files
- Run *cavity_prep.nb*

### Analysis

- Any of the notebooks can be run directly. 
- Tunable parameters are contained in the section **First run**, within subsections that are open by default
- During run-time, you will be prompted to enter additional parameters for *Algorithm 6.1* as more data becomes available.
  - Guidelines precede corresponding code segments.
- Notebooks are configured to save the output
  - Savefile can be used by commenting out the section **First run** and using **Later runs** 	




