# ceRvus
R package for interacting with the parentage analysis software [Cervus](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp) from within an R environment.

## Warning
You need to have an understanding of how Cervus works to use this package, as without the GUI, it is easy to misspecify an analysis and draw incorrect inferences from the results. This package is also only tested in a very narrow sense, with functions written primarily for use with my current project. It is made public to allow for others to replicate my results using the same scripts.

## Installation
The package can be currently installed from GitHub using the following commands:
```
# install.package("devtools")
devtools::install_github("irmoodie/ceRvus");
```
As Cervus is Windows only, the package will be of most use on Windows. To install the lastest version of Cervus (package developed using 3.0.7.0), follow the instructions on the [field genetics website](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp).

## Use
The package's goal is to reduce interaction with the Cervus GUI by providing a method to interact with the command line version of Cervus (CervusCL) from within an R environment. Currently, the package is being written to complement a meta-analysis I am working on, so development of features as and when I need them will take priority. This is also my first R package, so I am treating it as a learning exercise. I have made the code public in case someone finds it useful. **It is not tested beyond the needs of my current analysis.**

ceRvus can currently parameterise and run:

- Allele frequency analysis
- Parentage simulation
- Parentage analysis

from within R, as well as importing Cervus' output directly into R. It does this through creating/editing the .crv file that Cervus uses to store settings, and by parsing the .txt files that contain results into a format friendly for creating publication quality tables, or for performing further analysis within R.

In general, the package expects all the files for the analysis to be located within a single folder.

## Questions/Issues/Contributions
Open an [Issue](https://github.com/irmoodie/ceRvus/issues) on GitHub, and I will get back to you when I can.
