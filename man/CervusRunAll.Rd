% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CervusRunAll.R
\name{CervusRunAll}
\alias{CervusRunAll}
\title{Run Cervus ALF, SIM and PAR with an already configured project file}
\usage{
CervusRunAll(CervusCLPath, AnalysisFolderPath, ImportData = TRUE)
}
\arguments{
\item{CervusCLPath}{Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com).}

\item{AnalysisFolderPath}{Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.}

\item{ImportData}{TRUE/FALSE: Should the data be imported after the analysis has run?}
}
\value{
Outputs analysis summaries to the console, and if ImportData = TRUE returns a list with all data from the analysis.
}
\description{
Given an already configured Cervus project file, this function will re-run the full analysis. 
This is often helpful for writing fully reproducible analyses, or if you prefer to configure the analysis using the GUI of Cervus, then run from within R.
Ensure that only one Cervus project file is found within the AnalysisFolderPath (.crv file). I recommend having a folder within your R project dedicated to (each) Cervus analysis.
The function will update the paths for any files referenced in the cervus project file, if they do not match the current path of the analysis.
This function is NOT wine compatible
}
