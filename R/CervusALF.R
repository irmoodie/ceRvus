#' Perform an Allele Frequency Analysis using CervusCL
#'
#' Performs an Allele Frequency Analysis uisng CervusCL. 
#' You should be familiar with Cervus before using this function. 
#' The arguments relate directly to settings within Cervus, usually selected by the user using the GUI.
#' Argument names are such that it can easily be inferred what setting they change within a Cervus project file (.crv), often at the expense of brevity. This may change in future versions.
#' This function is wine compatible (see docs).
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus can be downloaded at http://www.fieldgenetics.com).
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder. I recommend making a folder within your project for this purpose.
#' @param AnalysisName User can specify a custom name for the analysis, which will be saved in the output filenames. Default is to use "CervusAnalysis".
#' @param ImportALF TRUE/FALSE: Import the alele frequency analysis summary tables in an R friendly format.
#' @param ResultsToConsole TRUE/FALSE: Should the results of the allele frequency analysis be output into the console?
#' @param GenotypeFile_FileName The filename of the genotype file. (e.g. genotype.csv).
#' @param GenotypeFile_HasHeader TRUE/FALSE: Does the genotype file contain a header row?
#' @param GenotypeFile_ReadLocusNames TRUE/FALSE: Should Cervus read the locus names in the file (TRUE), or give them generic names (FALSE).
#' @param GenotypeFile_IDColumnNumber Column number that identifies each individual.
#' @param GenotypeFile_FirstAlleleColumnNumber Column number which contains the first allele.
#' @param GenotypeFile_NLoci Number of loci in the genotype file.
#' @param DoHardyWeinberg TRUE/FALSE: Do Hardy-Weinberg test box to test each locus for whether it conforms to HW equilibrium.
#' @param HWMinExpectedFrequency Default is 5.
#' @param UseYatesCorrection TRUE/FALSE: Use Yates correction when calculating chi-square values with one degree of freedom.
#' @param UseBonferroniCorrection TRUE/FALSE: Use Bonferroni correction when carrying out the HW test for multiple loci simultaneously.
#' @param DoNullAllele TRUE/FALSE: Estimate null allele frequency.
#' @param wineCommand ONLY FOR WINE USERS: give the command that should be used to call your installed version of wine (e.g. wine64)
#' @param wineHomeDirectory ONLY FOR WINE USERS: Where should Wine find your home directory? This is usually "Z:".
#' @param wineTempDirectory ONLY FOR WINE USERS: Point wine to your Windows temp folder (e.g. "/Users/myname/.wine/drive_c/users/myname/Temp"). This is to fix a wine bug.
#' @param targetReturn Return output file paths for use with the targets package. Will break if ImportALF == TRUE
#' @return The function creates the .crv file that contains the settings for the analysis, then performs the analysis by supplying cervusCL.exe with the .crv file. The results are saved in the usual Cervus format, and printed to the console if desired. If ImportALF is set to TURE, the function also imports the summarised data into a list object. 
#' @examples 
#' CervusALF(
#'   CervusCLPath = "C:\\Program Files (x86)\\Field Genetics\\Cervus\\Cervus CL\\CervusCL.exe",
#'   AnalysisFolderPath = "C:\\Analysis",
#'   AnalysisName = "CervusAnalysis",
#'   ImportALF = TRUE,
#'   GenotypeFile_FileName = "genotype_file.csv",
#'   GenotypeFile_HasHeader = TRUE,
#'   GenotypeFile_ReadLocusNames = TRUE,
#'   GenotypeFile_IDColumnNumber = 1,
#'   GenotypeFile_FirstAlleleColumnNumber = 2,
#'   GenotypeFile_ColumnsPerLocus = 2,
#'   GenotypeFile_NLoci = 6,
#'   DoHardyWeinberg = TRUE,
#'   HWMinExpectedFrequency = 5,
#'   UseYatesCorrection = TRUE,
#'   UseBonferroniCorrection = TRUE,
#'   DoNullAllele = TRUE,
#'   ResultsToConsole = TRUE
# ')
#' @export

CervusALF <- 
  function(
    CervusCLPath,
    AnalysisFolderPath, 
    AnalysisName = "CervusAnalysis",
    ImportALF = FALSE,
    ResultsToConsole = FALSE,
    GenotypeFile_FileName,
    GenotypeFile_HasHeader = TRUE,
    GenotypeFile_ReadLocusNames = TRUE,
    GenotypeFile_IDColumnNumber = 1,
    GenotypeFile_FirstAlleleColumnNumber = 2,
    GenotypeFile_NLoci,
    DoHardyWeinberg = TRUE,
    HWMinExpectedFrequency = 5,
    UseYatesCorrection = TRUE,
    UseBonferroniCorrection = TRUE,
    DoNullAllele = TRUE,
    wineCommand = NA,
    wineHomeDirectory = "Z:",
    wineTempDirectory = NA,
    targetReturn = FALSE
  ) {
  
  # Specify paths to files
  pathGenotypeFile <- file.path(AnalysisFolderPath, GenotypeFile_FileName, fsep = .Platform$file.sep)
  pathAnalysisSettings <- file.path(AnalysisFolderPath, paste0(AnalysisName,"_settings", ".crv"), fsep = .Platform$file.sep)
  pathAlleleFrequencySummary <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_AlleleFrequencyAnalysis.txt"), fsep = .Platform$file.sep)
  pathAlleleFrequencyData <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_AlleleFrequencyAnalysis.alf"), fsep = .Platform$file.sep)
  
  # Check that CervusCLPath points to something
  if (!missing(CervusCLPath)) {
    if (file.exists(CervusCLPath)) {
      cervus_crv <- list(
        ProgramInfo = list(
          "ProgramName" = "Cervus",
          "ProgramVersion" = "3.0",
          "FileVersion" = "3.0.7.0",
          "CeRvusVersion" = paste0(packageVersion("ceRvus")),
          "Notes" = "This project file was generated with ceRvus, an R package to help interact with Cervus. Any issues related to it's usage should be reported at: https://github.com/irmoodie/ceRvus/issues"
        ),
        FileInfo = list(
          "FileName" = pathAnalysisSettings, # must be absolute path (e.g. C:\\...)
          "FileType" = ".crv",
          "CreationDate" = format(Sys.time(), usetz = FALSE)
        ),
        GenotypeFile = list(
          "FileName" = pathGenotypeFile, # must be absolute path (e.g. C:\\...)
          "HeaderRow" = paste(as.integer(GenotypeFile_HasHeader)), # If header row present, set to "1", if not set to "0"
          "ReadLocusNames" = paste(as.integer(GenotypeFile_ReadLocusNames)), # To read the names of loci, set to "1",
          "FirstAlleleColumnNumber" = paste0(GenotypeFile_FirstAlleleColumnNumber), # column number where alleles start,
          "IDColumnNumber" = paste0(GenotypeFile_IDColumnNumber), # column number that identifies each individual sampled,
          "NLoci" = paste0(GenotypeFile_NLoci), # number of loci to be used
          "SexColumn" = "0", # should probably remain as "0", don't think it would be required here
          "UnknownSexLabel" = "" # not required at this stage
        ),
        AlleleFrequencySummaryFile = list(
          "FileName" = pathAlleleFrequencySummary, # must be absolute path (e.g. C:\\...)
          "DoHardyWeinberg" = paste(as.integer(DoHardyWeinberg)), # "1" to do Hardy-Weinberg test box to test each locus for whether it conforms to HW equilibrium
          "HWMinExpectedFrequency" = paste(as.integer(HWMinExpectedFrequency)), # "5" is default
          "UseYatesCorrection" = paste(as.integer(UseYatesCorrection)), # "1" to use Yates correction when calculating chi-square values with one degree of freedom
          "UseBonferroniCorrection" = paste(as.integer(UseBonferroniCorrection)), # "1" to use Bonferroni correction when carrying out the HW test for multiple loci simultaneously
          "DoNullAllele" = paste(as.integer(DoNullAllele)) # "1" to estimate null allele frequency
        ),
        AlleleFrequencyDataFile = list(
          "FileName" = pathAlleleFrequencyData, # must be absolute path (e.g. C:\\...)
          "HeaderRow" = "1"
        ))
      
      if (!is.na(wineCommand)) {
        cervus_crv$FileInfo$FileName <- paste0(wineHomeDirectory, pathAnalysisSettings)
        cervus_crv$GenotypeFile$FileName = paste0(wineHomeDirectory, pathGenotypeFile)
        cervus_crv$AlleleFrequencySummaryFile$FileName = paste0(wineHomeDirectory, pathAlleleFrequencySummary)
        cervus_crv$AlleleFrequencyDataFile$FileName = paste0(wineHomeDirectory, pathAlleleFrequencyData)
        
        # Normal cervus runs perfect, but currently a wine bug is breaking the command line version
        # https://bugs.winehq.org/show_bug.cgi?id=49334
        # as a work around, r will create the temp file locations needed
        # you shouldn't have to worry about temp folder bloat, as CervusCL clears it after it runs
        
        dir.create(
          path = file.path(wineTempDirectory, AnalysisFolderPath, fsep = .Platform$file.sep), 
          recursive = TRUE, 
          showWarnings = FALSE)
        
      }
      
      
      ini::write.ini(x = cervus_crv, filepath = pathAnalysisSettings) # requires ini package to format settings file
      
      # No Wine run
      
      if (is.na(wineCommand)) {
        system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/ALF /O")) # run the allele frequency analysis
        if (ResultsToConsole) {
          system(command = paste0("cat ", '"', pathAlleleFrequencySummary, '"')) # display the results 
        }
      }
      
      # Wine run
      
      if (!is.na(wineCommand)) {
        system(command = paste0(wineCommand, ' "', CervusCLPath, '" ', '"', wineHomeDirectory, pathAnalysisSettings, '" ', "/ALF /O")) # run the allele frequency analysis
        if (ResultsToConsole) {
        system(command = paste0("cat ", '"', pathAlleleFrequencySummary, '"'))  # display the results
        }
      }

      if (ImportALF) {
        ALFSummary <- ImportCervusALF(ALFSummaryFile = pathAlleleFrequencySummary)
        return(ALFSummary)
      }
      
      cat("\nAnalysis complete!\n")
      
      if (targetReturn) {
        output_files <- c(pathAnalysisSettings, pathAlleleFrequencySummary, pathAlleleFrequencyData)
        return(output_files)
      }
      
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  } else {
    cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
  }
  
}
