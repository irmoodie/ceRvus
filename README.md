<p align="center">
<img src="https://user-images.githubusercontent.com/38960700/211855393-96d482f1-e0ba-4393-bc83-328a69497111.png" width="259" height="300">
</p>

<p align="center" href="https://doi.org/10.5281/zenodo.17968976"><img src="https://zenodo.org/badge/498382224.svg" alt="DOI"></p>

# ceRvus
R package for interacting with the parentage analysis software [Cervus](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp) from within an R environment. The package also includes functions to make working with CervusCL via Wine (i.e. on Linux/MacOS) easier.

## Warning
You need to have an understanding of how Cervus works to use this package fully, as without the GUI, it is easy to misspecify an analysis and draw incorrect inferences from the results. However, it is also possible to run an analysis that has been specified using the GUI of Cervus, which may be preferable depending on your workflow. This package is also only tested in a narrow sense, with functions written primarily for use with my current project. It is made public to allow for others to replicate my results using the same scripts, and to help others who might want to do similar things.

## Installation
The package can be currently installed from GitHub using the following commands:
```
# install.package("devtools")
devtools::install_github("irmoodie/ceRvus");
```
To install the lastest version of Cervus (package developed using 3.0.7.0), follow the instructions on the [field genetics website](http://www.fieldgenetics.com/pages/aboutCervus_Overview.jsp).

The package relies on the [ini](https://github.com/dvdscripter/ini), [dplyr](https://github.com/tidyverse/dplyr), [stringr](https://github.com/tidyverse/stringr) and [tidyr](https://github.com/tidyverse/tidyr) packages.

### Using MacOS or Linux?

It is quite easy to run Cervus via [Wine](https://www.winehq.org/) on both MacOS and Linux. I've tested that the results are identical across platforms, and find run times are as expected (and much faster compared to running a VM). For more info see: [Running Cervus using Wine](https://github.com/irmoodie/ceRvus/blob/master/docs/CervusWine.md)

## Use
The package's goal is to reduce interaction with the Cervus GUI by providing a method to interact with the command line version of Cervus (CervusCL) from within an R environment. Currently, the package is being written to complement a meta-analysis I am working on, so development of features as and when I need them will take priority. This is also my first R package, so I am treating it as a learning exercise. I have made the code public in case someone finds it useful. **It is not tested beyond the needs of my current analysis.**

ceRvus can currently parameterise and run:

- Allele frequency analysis
- Parentage simulation
- Parentage analysis

(It currently cannot handle identity analyses, but this would be easy to add if required)

from within R, as well as importing Cervus' output directly into R. It does this through creating/editing the .crv file that Cervus uses to store settings, and by parsing the .txt files that contain results into a format friendly for creating publication quality tables, or for performing further analysis within R.

In general, the package expects all the files for the analysis to be located within a single folder.

It is left up to the user to decide if they want to work with their genotype data in R, then export, or just work with already existing datafiles. Care should be taken when exporting from R to then be used in Cervus, ESPECIALLY WITH REGARD TO NAs. Cervus will treat the string 'NA' as an allele and estimate it's frequency etc!

The package is by no-means perfect, it is very janky and hacky, but it functions. As all the heavy lifting is still done by CervusCL, so there is no performance loss or gain. 

## Citation
If you find ceRvus helpful, you can cite the following DOI: [https://doi.org/10.5281/zenodo.17968976](https://doi.org/10.5281/zenodo.17968976)

Or use the following BibLaTeX entry:

```tex
@software{moodieCervus2026,
  title = {ceRvus v0.0.3.0},
  author = {Moodie, Iain R.},
  date = {2026-01-13},
  doi = {10.5281/zenodo.18231860},
  url = {https://github.com/irmoodie/ceRvus},
}
```

## Questions/Issues/Contributions
Open an [Issue](https://github.com/irmoodie/ceRvus/issues) on GitHub, and I will get back to you when I can. Any suggestions, please open a pull request. 
