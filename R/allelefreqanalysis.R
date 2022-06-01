#' Allele Frequency Analysis
#'
#' Perform an Allele Frequency Analysis uisng CervusCL
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com)
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param AnalysisName Can specify a custom name for the analysis, which will be saved as the .crv filename
#' @param GenotypeFile_FileName The filename of the genotype file. (e.g. genotype.csv)
#' @param GenotypeFile_HasHeader TRUE/FALSE: Does the genotype file contain a header row
#' @param GenotypeFile_ReadLocusNames TRUE/FALSE: Should Cervus read the locus names in the file (TRUE), or give them generic names (FALSE)
#' @param GenotypeFile_IDColumnNumber Column number that identifies each individual
#' @param GenotypeFile_FirstAlleleColumnNumber Column number which contains the first allele
#' @param GenotypeFile_ColumnsPerLocus What format are the loci in? One allele per column = 2 columns per loci. See Cervus helpfiles for more info.
#' @param GenotypeFile_NLoci Number of loci in the genotype file
#' @param DoHardyWeinberg TRUE/FALSE: Do Hardy-Weinberg test box to test each locus for whether it conforms to HW equilibrium
#' @param HWMinExpectedFrequency Default is 5
#' @param UseYatesCorrection TRUE/FALSE: Use Yates correction when calculating chi-square values with one degree of freedom
#' @param UseBonferroniCorrection TRUE/FALSE: Use Bonferroni correction when carrying out the HW test for multiple loci simultaneously
#' @param DoNullAllele TRUE/FALSE Estimate null allele frequency
#' @return Currently the function creates the .crv file that contains the settings for the analysis, then performs the analysis. The results are saved in the usual CERVUS format, and printed to the console.
#' @examples 
#' CervusALF(
#'   CervusCLPath = "C:\\Program Files (x86)\\Field Genetics\\Cervus\\Cervus CL\\CervusCL.exe",
#'   AnalysisFolderPath = "C:\\Analysis",
#'   AnalysisName = "CERVUS_Analysis",
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
#'   DoNullAllele = TRUE
# ')
#' @export

CervusALF <- function(CervusCLPath,
                       AnalysisFolderPath, 
                       AnalysisName = "CERVUS_Analysis", 
                       GenotypeFile_FileName,
                       GenotypeFile_HasHeader = TRUE,
                       GenotypeFile_ReadLocusNames = TRUE,
                       GenotypeFile_IDColumnNumber = 1,
                       GenotypeFile_FirstAlleleColumnNumber = 2,
                       GenotypeFile_ColumnsPerLocus = 2,
                       GenotypeFile_NLoci,
                       DoHardyWeinberg = TRUE,
                       HWMinExpectedFrequency = 5,
                       UseYatesCorrection = TRUE,
                       UseBonferroniCorrection = TRUE,
                       DoNullAllele = TRUE
) {
  
  # Specify paths to files
  pathGenotypeFile <- file.path(AnalysisFolderPath, GenotypeFile_FileName, fsep = "\\")
  pathAnalysisSettings <- file.path(AnalysisFolderPath, paste0(AnalysisName,"_settings", ".crv"), fsep = "\\")
  pathAlleleFrequencySummary <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_AlleleFrequencyAnalysis.txt"), fsep = "\\")
  pathAlleleFrequencyData <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_AlleleFrequencyAnalysis.alf"), fsep = "\\")
  
  # Check that CervusCLPath points to something
  if (!missing(CervusCLPath)) {
    if (file.exists(CervusCLPath)) {
      cervus_crv <- list(
        ProgramInfo = list(
          "ProgramName" = "Cervus",
          "ProgramVersion" = "3.0",
          "FileVersion" = "3.0.7.0"
        ),
        Registration = list(
          "UserName" = "",
          "UserCompany" = "",
          "Code" = ""
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
          "PropLociTyped" = "0", # should remain as "0" here,
          "ColumnsPerLocus" = paste0(GenotypeFile_ColumnsPerLocus), # set as "2" assuming the default input format for CERVUS
          "SexColumn" = "0", # should probably remain as "0", don't think it would be required here
          "UnknownSexLabel" = ""
        ),
        CodecFile = list(
          "FileName" = "", # must be absolute path (e.g. C:\\...)
          "HeaderRow" = "1", # If header row present, set to "1", if not set to "0"
          "UseSameCodingForAllLoci" = "1", # should probably be kept as "1"
          "GenotypeFileName" = pathGenotypeFile # must be absolute path (e.g. C:\\...)
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
          "HeaderRow" = "1" # If header row present, set to "1", if not set to "0"
        ))
      
      ini::write.ini(x = cervus_crv, filepath = cervus_crv$FileInfo$FileName) # requires ini package to format settings file
      
      system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/ALF /O")) # run the allele frequency analysis
      system(command = paste0("cat ", '"', pathAlleleFrequencySummary, '"')) # display the results
      
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  } else {
    cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
  }
  
}