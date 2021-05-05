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
- Tunable parameters are contained in **Computations**/**Temporal parameters**, within subsections that are open by default
- During run-time, you will be prompted to enter additional parameters for *Algorithm 6.1* as more data becomes available.
  - Guidelines precede corresponding code segments.
- Notebooks are configured to save the output, and hence restart the analysis, at two distinct stages
  - **Restore-point #1**: After DMD over the chosen set of delays, model orders and trajectories, but before any further post-processing
  - **Restore-point #2**: At the end of all intended analysis
- By default, the code saves output at both stages so that one can step back to either of those points
- Restore procedure: Run **Initialization**, then **Restore-point** of choice, and proceed as usual. 







