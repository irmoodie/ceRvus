#' Run Cervus ALF, SIM and PAR with an already configured project file
#'
#' Given an already configured Cervus project file, this function will re-run the full analysis. 
#' This is often helpful for writing fully reproducible analyses, or if you prefer to configure the analysis using the GUI of Cervus, then run from within R.
#' Ensure that only one Cervus project file is found within the AnalysisFolderPath (.crv file). I recommend having a folder within your R project dedicated to (each) Cervus analysis.
#' The function will update the paths for any files referenced in the cervus project file, if they do not match the current path of the analysis.
#' This function is NOT wine compatible
#' 
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com).
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param ImportData TRUE/FALSE: Should the data be imported after the analysis has run?
#' @return Outputs analysis summaries to the console, and if ImportData = TRUE returns a list with all data from the analysis.
#' @export

CervusRunAll <- function(CervusCLPath, AnalysisFolderPath, ImportData = TRUE){
  
  if (!missing(CervusCLPath)) {
    if (file.exists(CervusCLPath)) {
      
      pathAnalysisSettings <- Sys.glob(file.path(AnalysisFolderPath, "*.crv", fsep = .Platform$file.sep))
      
      cervus_crv_file <- ini::read.ini(pathAnalysisSettings)
      
      change_detect <- 0
      
      for (var in c("FileInfo", "GenotypeFile", "CodecFile", "AlleleFrequencySummaryFile", "AlleleFrequencyDataFile", "SimulationSummaryFile", "SimulationDataFile", "OffspringFile", "CandidateFemaleFile", "CandidateMaleFile", "ParentageSummaryFile", "ParentageDataFile")) {
        if (var != "CodecFile") {
          if (exists(x = "FileName", where = cervus_crv_file[[var]])){
            if (cervus_crv_file[[var]]$FileName != file.path(AnalysisFolderPath, basename(cervus_crv_file[[var]]$FileName),  fsep = .Platform$file.sep)) {
              cervus_crv_file[[var]]$FileName <- file.path(AnalysisFolderPath, basename(cervus_crv_file[[var]]$FileName),  fsep = .Platform$file.sep)
              change_detect <- change_detect+1
            }
          }
        }
        if (var == "CodecFile") {
          if (exists(x = "GenotypeFileName", where = cervus_crv_file[[var]])){
            if (cervus_crv_file$CodecFile$GenotypeFileName != file.path(AnalysisFolderPath, basename(cervus_crv_file$CodecFile$GenotypeFileName),  fsep = .Platform$file.sep)) {
              cervus_crv_file$CodecFile$GenotypeFileName <- file.path(AnalysisFolderPath, basename(cervus_crv_file$CodecFile$GenotypeFileName),  fsep = .Platform$file.sep)
              change_detect <- change_detect+1
            }
          }
        }
      }
      
      if (change_detect != 0) {
        ini::write.ini(x = cervus_crv_file, filepath = pathAnalysisSettings)
      }
      
      system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/ALF /O")) # run the allele frequency analysis
      system(command = paste0("cat ", '"', cervus_crv_file$AlleleFrequencySummaryFile$FileName, '"')) # display the results
      
      system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/SIM /O")) # run the simulation
      system(command = paste0("cat ", '"', cervus_crv_file$SimulationSummaryFile$FileName, '"')) # display the results
      
      system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/PAR /O")) # run the analysis
      system(command = paste0("cat ", '"', cervus_crv_file$ParentageSummaryFile$FileName, '"')) # display the results
      
      if (isTRUE(ImportData)) {
        data <- list()
        data$AlleleFrequencyAnalysis <- ImportCervusALF(cervus_crv_file$AlleleFrequencySummaryFile$FileName)
        return(data)
      }
      
      cat("\nAnalysis complete!")
      
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  } else {
    cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
  }
}
