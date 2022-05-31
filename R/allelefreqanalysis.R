#' Allele Frequency Analysis
#'
#' Perform an Allele Frequency Analysis uisng CervusCL
#' @import ini
#' @param cervusCL_location Path to cervusCL.exe (CERVUS is found at http://www.fieldgenetics.com)
#' @param analysis_folder Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param analysis_name Can specify a custom name for the analysis, which will be saved as the .crv filename
#' @param genotype_file The filename of the genotype file. Currently, the function expects a .csv file with the first column as the id, then the loci starting in column two. Headers are assumed and loci names are retained.
#' @param NLoci Number of loci in the genotype file
#' @param DoHardyWeinberg TRUE/FALSE Do Hardy-Weinberg test box to test each locus for whether it conforms to HW equilibrium
#' @param HWMinExpectedFrequency Default is 5
#' @param UseYatesCorrection TRUE/FALSE Use Yates correction when calculating chi-square values with one degree of freedom
#' @param UseBonferroniCorrection TRUE/FALSE Use Bonferroni correction when carrying out the HW test for multiple loci simultaneously
#' @param DoNullAllele TRUE/FALSE Estimate null allele frequency
#' @return Currently the function creates the .crv file that contains the settings for the analysis, then performs the analysis. The results are saved in the usual CERVUS format, and printed to the console.
#' @examples 
#' cervus_alf(cervusCL_location = cervusCL_location,
#'            analysis_folder = analysis_folder,
#'            genotype_file = "genotype_file.csv",
#'            analysis_name = "test_analysis",
#'            NLoci = 6,
#'            DoHardyWeinberg = TRUE,
#'            HWMinExpectedFrequency = 5,
#'            UseYatesCorrection = TRUE,
#'            UseBonferroniCorrection = TRUE,
#'            DoNullAllele = TRUE
#'            )
#' 
#' @export

cervus_alf <- function(cervusCL_location,
                       analysis_folder, 
                       analysis_name = "CERVUS_Analysis", 
                       genotype_file,
                       NLoci,
                       DoHardyWeinberg = TRUE,
                       HWMinExpectedFrequency = 5,
                       UseYatesCorrection = TRUE,
                       UseBonferroniCorrection = TRUE,
                       DoNullAllele = TRUE
) {
  
  GenotypeFileName <- file.path(analysis_folder, genotype_file, fsep = "\\")
  AnalysisSettingsCrv <- file.path(analysis_folder, paste0(analysis_name, ".crv"), fsep = "\\")
  AlleleFrequencySummaryFileName <- file.path(analysis_folder, "AlleleFrequencyAnalysis.txt", fsep = "\\")
  AlleleFrequencyDataFileName <- file.path(analysis_folder, "AlleleFrequencyAnalysis.alf", fsep = "\\")
  
  if (!missing(cervusCL_location)) {
    if (file.exists(cervusCL_location)) {
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
          "FileName" = AnalysisSettingsCrv, # must be absolute path (e.g. C:\\...)
          "FileType" = ".crv",
          "CreationDate" = format(Sys.time(), usetz = FALSE)
        ),
        GenotypeFile = list(
          "FileName" = GenotypeFileName, # must be absolute path (e.g. C:\\...)
          "HeaderRow" = "1", # If header row present, set to "1", if not set to "0"
          "ReadLocusNames" = "1", # To read the names of loci, set to "1",
          "FirstAlleleColumnNumber" = "2", # column number where alleles start,
          "IDColumnNumber" = "1", # column number that identifies each individual sampled,
          "NLoci" = paste0(NLoci), # number of loci to be used
          "PropLociTyped" = "0", # should remain as "0" here,
          "ColumnsPerLocus" = "2", # set as "2" assuming the default input format for CERVUS
          "SexColumn" = "0", # should probably remain as "0", don't think it would be required here
          "UnknownSexLabel" = ""
        ),
        CodecFile = list(
          "FileName" = "", # must be absolute path (e.g. C:\\...)
          "HeaderRow" = "1", # If header row present, set to "1", if not set to "0"
          "UseSameCodingForAllLoci" = "1", # should probably be kept as "1"
          "GenotypeFileName" = GenotypeFileName # must be absolute path (e.g. C:\\...)
        ),
        AlleleFrequencySummaryFile = list(
          "FileName" = AlleleFrequencySummaryFileName, # must be absolute path (e.g. C:\\...)
          "DoHardyWeinberg" = paste(as.integer(DoHardyWeinberg)), # "1" to do Hardy-Weinberg test box to test each locus for whether it conforms to HW equilibrium
          "HWMinExpectedFrequency" = paste(as.integer(HWMinExpectedFrequency)), # "5" is default
          "UseYatesCorrection" = paste(as.integer(UseYatesCorrection)), # "1" to use Yates correction when calculating chi-square values with one degree of freedom
          "UseBonferroniCorrection" = paste(as.integer(UseBonferroniCorrection)), # "1" to use Bonferroni correction when carrying out the HW test for multiple loci simultaneously
          "DoNullAllele" = paste(as.integer(DoNullAllele)) # "1" to estimate null allele frequency
        ),
        AlleleFrequencyDataFile = list(
          "FileName" = AlleleFrequencyDataFileName, # must be absolute path (e.g. C:\\...)
          "HeaderRow" = "1" # If header row present, set to "1", if not set to "0"
        ))
      
      if (require(ini)) {
        ini::write.ini(x = cervus_crv, filepath = cervus_crv$FileInfo$FileName) # requires ini package to format settings file
      } else {
        cat("WARNING: Please install the 'ini' package using install.package('ini')")
      }
      
      system(command = paste0('"', cervusCL_location, '" ', '"', AnalysisSettingsCrv, '" ', "/ALF /O")) # run the allele frequency analysis
      system(command = paste0("cat ", '"', AlleleFrequencySummaryFileName, '"')) # display the results
      
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  } else {
    cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
  }
  
}