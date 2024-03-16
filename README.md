# Dynamical-Tidal-Response-of-Polytropic-Stars

This repository contains a code to analyze the dynamical tidal response of polytropic stars in general relativity.
The code is based on the paper [Dynamical tidal response of non-rotating relativistic stars](https://arxiv.org/abs/2403.03254). If you use this code please use the citation
```
@article{HegadeKR:2024agt,
    author = "Hegade K. R., Abhishek and Ripley, Justin L. and Yunes, Nicol\'as",
    title = "{Dynamical tidal response of non-rotating relativistic stars}",
    eprint = "2403.03254",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "3",
    year = "2024"
}
```

## Description of the Mathematica Notebooks

The repository contains a number of Mathematica notebooks.

The notebook supplementary_funcs.nb contains the expressions for the supplementary functions defined in the paper.

The notebook generate-normalized-solution.nb generates the normalized solutions $\hat{N} H_{P/Q,\ell}$ used in the paper for a given value of $\ell$.

## Description of C++ Code
The C++ code provided in the folder Polytropic-EoS can be used to study the dynamical tidal response of a polytropic star in full general relativity.
To install the code you need to install the C++ package [boost](https://formulae.brew.sh/formula/boost) and download the eigen-3.4.0 package from [here](https://eigen.tuxfamily.org/index.php?title=Main_Page).

After installing the above packages please change the BOOST_ROOT line and the EIGEN_ROOT line in the Makefile. 
Once this is set the code should be ready to install. 
First create a folder called obj 
```
mkdir obj
```
then use 
```
make
```
the make file should compile without any errors.

To run the code type 
```
setup_run.py 
```

The output will be saved in the Runs folder. 
