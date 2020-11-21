# DIMER_INTERACTIONS
FF vs QM comparison for interactions between 2 molecules

# Requirements:
Install anaconda or miniconda and run the following:
```
  conda create -n gmx_plumed python=2.7  #(some of the older alchemical_analysis packages work only with 2.7)
  conda activate gmx_plumed  
  conda install -c ajaymur91 -c conda-forge plumed   # plumed 2.6.1
  conda install -c ajaymur91 -c conda-forge gromacs  # gromacs 2019.6, plumed patched version, using conda-build and passed 37/37 regression tests, simd=sse4.1, no mpi
  conda install -c openbabel openbabel
  ```
# Usage
  ```
  conda activate gmx_plumed 
  git clone https://github.com/ajaymur91/DIMER_INTERACTIONS.git
  cd DIMER_INTERACTIONS
  bash run_n.sh
  ```
  - Creates dimer configurations from gas phase MD -> gro/
  - Saves FF interaction energies of those configurations in results/E_FF_kcal.txt
  - Converts gro files to gaussian input files -> com/
  - Currently supports LI, TFSI, BLB and EC only

