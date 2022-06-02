#' Run Cervus ALF, SIM and PAR with an already configured project file
#'
#' Given an already configured Cervus project file, this function will re-run the full analysis. 
#' This is often helpful for writing fully reproducible analyses.
#' Ensure that only one Cervus project file is found within the AnalysisFolderPath (.crv file).
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com).
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param ImportData TRUE/FALSE: Should the data be imported after the analysis has run?
#' @return Outputs analysis summaries to the console, and if ImportData = TRUE returns a list with all data from the analysis.
#' @export

CervusRunAll <- function(CervusCLPath, AnalysisFolderPath, ImportData = TRUE){
  
  pathAnalysisSettings <- Sys.glob(file.path(AnalysisFolderPath, "*.crv", fsep = "\\"))
  
  CervusCRVFile <- ini::read.ini(pathAnalysisSettings)
  
  system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/ALF /O")) # run the allele frequency analysis
  system(command = paste0("cat ", '"', CervusCRVFile$AlleleFrequencySummaryFile$FileName, '"')) # display the results
  
  system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/SIM /O")) # run the simulation
  system(command = paste0("cat ", '"', CervusCRVFile$SimulationSummaryFile$FileName, '"')) # display the results
  
  system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/PAR /O")) # run the analysis
  system(command = paste0("cat ", '"', CervusCRVFile$ParentageSummaryFile$FileName, '"')) # display the results
  
  cat("\nAnalysis complete!")
  
  if (isTRUE(ImportData)) {
    data <- list()
    data$AlleleFrequencyAnalysis <- ImportCervusALF(CervusCRVFile$AlleleFrequencySummaryFile$FileName)
    return(data)
    
    cat("\nData import complete!")
  }
  
}
