# ceRvus
An R package for interacting with the parentage analysis software [Cervus](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp) from within an R environment.

## Installation
The package can be currently installed from GitHub using the following commands:
```
# install.package("devtools")
devtools::install_github("irmoodie/ceRvus");
```
As Cervus is Windows only, the package will be of most use on Windows. Future updates may include functions that would be useful to users who have the output of a Cervus analysis that they want to work with in R. To install the lastest version of Cervus (package developed using 3.0.7.0), follow the instructions on the [field genetics website](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp).

## Use
The package's goal is to reduce interaction with the Cervus GUI by providing a method to interact with the command line version of Cervus (CervusCL) from within an R environment.
Currently, the package is being written to complement a meta-analysis I am working on, so development of features as and when I need them will take priority.
This is also my first R package, so I am treating it as a learning exercise. I have made the code public in case someone finds it useful.

Current goals for the package are to be able to parameterise and run:

- Allele frequency analysis
- Parentage simulation
- Parentage analysis

from within R, as well as importing Cervus' output directly into R.
