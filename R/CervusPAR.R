#' Run parentage analysis using CervusCL
#'
#' You should be familiar with Cervus before using this function. 
#' The arguments relate directly to settings within Cervus, usually selected by the user using the GUI.
#' Argument names are such that it can easily be inferred what setting they change within a Cervus project file (.crv), often at the expense of brevity. This may change in future versions.
#' This function is wine compatible (see docs).
#' 
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com).
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param AnalysisName Can specify a custom name for the analysis, which will be saved in the output filenames.
#' @param AnalysisType Set as "Maternity"  or "Paternity" or "Parent pair (sexes known)" or "Parent pair (sexes unknown)".
#' @param OffspringFile_FileName Name of the file with offspring data
#' @param OffspringFile_HasHeader TRUE/FALSE: Does the file have a header?
#' @param OffspringFile_OffspringIDColumnNumber Column number with offspring ID
#' @param OffspringFile_IncludesKnownParents TRUE/FALSE
#' @param OffspringFile_KnownParentIDColumnNumber Column number with offspring ID (set to 0 if no known parents)
#' @param OffspringFile_IncludesCandidateParents TRUE/FALSE: does the offspring file also contain candidate parents?
#' @param OffspringFile_CandidateParentIDColumnNumber Column number where candidate parent IDs start (set to 0 if no candidate parents in the file)
#' @param CandidateFemaleFile_FileName Name of the file with candidate females
#' @param CandidateFemaleFile_HasHeader TRUE/FALSE: Does the file have a header?
#' @param CandidateFemaleFile_CandidateParentFormat "One column for all offspring"
#' @param CandidateFemaleFile_OffspringIDColumnNumber If using a format other than "One column for all offspring", give the column number that contains the offspring ID
#' @param CandidateFemaleFile_CandidateParentIDColumnNumber If using a format other than "One column for all offspring", give the column number that contains the first candidate parent ID
#' @param CandidateMaleFile_FileName Name of the file with candidate males
#' @param CandidateMaleFile_HasHeader TRUE/FALSE: Does the file have a header?
#' @param CandidateMaleFile_CandidateParentFormat "One column for all offspring"
#' @param CandidateMaleFile_OffspringIDColumnNumber If using a format other than "One column for all offspring", give the column number that contains the offspring ID
#' @param CandidateMaleFile_CandidateParentIDColumnNumber If using a format other than "One column for all offspring", give the column number that contains the first candidate parent ID
#' @param UseSimulationParameters TRUE/FALSE: TRUE (default) to use the SIM parameters
#' @param CalculateConfidenceLevels TRUE/FALSE: TRUE (default)
#' @param AlwaysTestSelfing TRUE/FALSE: Should selfing be tested?
#' @param MinTypedLoci Minimum number of typed loci to be included. Should match that of the simulation.
#' @param UseCorrectedLikelihoods TRUE/FALSE: Use Kalinowski et al. (2007) likelihood equations (recommended)
#' @param UseMistypingRateAsLikelihoodErrorRate TRUE/FALSE
#' @param LikelihoodErrorRate Default is 0.01
#' @param CriticalStatisticName Delta" or "LOD". "Delta" recommended. See the Cervus help files.
#' @param TruncateAtZero TRUE/FALSE
#' @param OutputType Default is "The most-likely parent", see the Cervus help files for more.
#' @param OutputSortedBy Default is "Joint LOD score", see the Cervus help files for more.
#' @param OutputIncludeNonExclusionProbabilities TRUE/FALSE
#' @param wineCommand ONLY FOR WINE USERS: give the command that should be used to call your installed version of wine (e.g. wine64)
#' @param wineHomeDirectory ONLY FOR WINE USERS: Where should Wine find your home directory? This is usually "Z:".
#' @param wineTempDirectory ONLY FOR WINE USERS: Point wine to your Windows temp folder (e.g. "/Users/myname/.wine/drive_c/users/myname/Temp"). This is to fix a wine bug.
#' @param ResultsToConsole TRUE/FALSE: Should the results of the allele frequency analysis be output into the console?
#' @return Currently prints the results to the console, and saves the parentage analysis data in AnalysisFolderPath
#'
#' @export

CervusPAR <- 
  function(
    CervusCLPath,
    AnalysisFolderPath, 
    AnalysisName = "CervusAnalysis",
    AnalysisType,
    OffspringFile_FileName,
    OffspringFile_HasHeader = TRUE,
    OffspringFile_OffspringIDColumnNumber = 1,
    OffspringFile_IncludesKnownParents = FALSE,
    OffspringFile_KnownParentIDColumnNumber = 0,
    OffspringFile_IncludesCandidateParents = FALSE,
    OffspringFile_CandidateParentIDColumnNumber = 0,
    CandidateFemaleFile_FileName = "",
    CandidateFemaleFile_HasHeader = FALSE,
    CandidateFemaleFile_CandidateParentFormat = "One column for all offspring",
    CandidateFemaleFile_OffspringIDColumnNumber = "",
    CandidateFemaleFile_CandidateParentIDColumnNumber = "",
    CandidateMaleFile_FileName = "",
    CandidateMaleFile_HasHeader = FALSE,
    CandidateMaleFile_CandidateParentFormat = "One column for all offspring",
    CandidateMaleFile_OffspringIDColumnNumber = "",
    CandidateMaleFile_CandidateParentIDColumnNumber = "",
    UseSimulationParameters = TRUE,
    CalculateConfidenceLevels = TRUE,
    AlwaysTestSelfing = FALSE,
    MinTypedLoci,
    UseCorrectedLikelihoods = TRUE,
    UseMistypingRateAsLikelihoodErrorRate = TRUE,
    LikelihoodErrorRate = 0.01,
    CriticalStatisticName = "Delta",
    TruncateAtZero = TRUE,
    OutputType = "The most-likely parent",
    OutputSortedBy = "Joint LOD score",
    OutputIncludeNonExclusionProbabilities = FALSE,
    ResultsToConsole = TRUE,
    wineCommand = NA,
    wineHomeDirectory = "Z:",
    wineTempDirectory = NA
  ){
    
    if (!missing(CervusCLPath)) {
      if (file.exists(CervusCLPath)) {
        
        pathAnalysisSettings <- file.path(AnalysisFolderPath, paste0(AnalysisName,"_settings", ".crv"), fsep = .Platform$file.sep)
        pathParentageSummaryFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Parentage.txt"), fsep = .Platform$file.sep)
        pathParentageDataFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Parentage.sim"), fsep = .Platform$file.sep)
        
        pathOffspringFile <- file.path(AnalysisFolderPath, OffspringFile_FileName, fsep = .Platform$file.sep)
        
        if (CandidateFemaleFile_FileName != "") {
          pathCandidateFemaleFile <- file.path(AnalysisFolderPath, CandidateFemaleFile_FileName, fsep = .Platform$file.sep)
        } else {
          pathCandidateFemaleFile <- CandidateFemaleFile_FileName
        }
        
        if (CandidateMaleFile_FileName != "") {
          pathCandidateMaleFile <- file.path(AnalysisFolderPath, CandidateMaleFile_FileName, fsep = .Platform$file.sep)
        } else {
          pathCandidateMaleFile <- CandidateMaleFile_FileName
        }
        
        cervus_crv_file <- ini::read.ini(pathAnalysisSettings)
        
        cervus_crv_par <- list(
          ParentageParameters = list(
            "AnalysisType" = paste0(AnalysisType),
            "UseSimulationParameters" = paste0(as.integer(UseSimulationParameters)),
            "CalculateConfidenceLevels" = paste0(as.integer(CalculateConfidenceLevels)),
            "AlwaysTestSelfing" = paste0(as.integer(AlwaysTestSelfing)),
            "MinTypedLoci" = paste0(MinTypedLoci),
            "UseCorrectedLikelihoods" = paste0(as.integer(UseCorrectedLikelihoods)),
            "UseMistypingRateAsLikelihoodErrorRate" = paste0(as.integer(UseMistypingRateAsLikelihoodErrorRate)),
            "LikelihoodErrorRate" = paste0(LikelihoodErrorRate),
            "CriticalStatisticName" = paste0(CriticalStatisticName),
            "TruncateAtZero" = paste0(as.integer(TruncateAtZero))
          ),
          OffspringFile = list(
            "FileName" = pathOffspringFile,
            "HeaderRow" = paste0(as.integer(OffspringFile_HasHeader)),
            "OffspringIDColumnNumber" = paste0(OffspringFile_OffspringIDColumnNumber),
            "IncludesKnownParents" = paste0(as.integer(OffspringFile_IncludesKnownParents)),
            "KnownParentIDColumnNumber" = paste0(OffspringFile_KnownParentIDColumnNumber),
            "IncludesCandidateParents" = paste0(as.integer(OffspringFile_IncludesCandidateParents)),
            "CandidateParentIDColumnNumber" = paste0(OffspringFile_CandidateParentIDColumnNumber)
          ),
          CandidateFemaleFile = list(
            "FileName" = pathCandidateFemaleFile,
            "HeaderRow" = paste0(as.integer(CandidateFemaleFile_HasHeader)),
            "CandidateParentFormat" = paste0(CandidateFemaleFile_CandidateParentFormat),
            "OffspringIDColumnNumber" = paste0(CandidateFemaleFile_OffspringIDColumnNumber),
            "CandidateParentIDColumnNumber" = paste0(CandidateFemaleFile_CandidateParentIDColumnNumber)
          ),
          CandidateMaleFile = list(
            "FileName" = pathCandidateMaleFile,
            "HeaderRow" = paste0(as.integer(CandidateMaleFile_HasHeader)),
            "CandidateParentFormat" = paste0(CandidateMaleFile_CandidateParentFormat),
            "OffspringIDColumnNumber" = paste0(CandidateMaleFile_OffspringIDColumnNumber),
            "CandidateParentIDColumnNumber" = paste0(CandidateMaleFile_CandidateParentIDColumnNumber)
          ),
          ParentageSummaryFile = list(
            "FileName" = pathParentageSummaryFile
          ),
          ParentageDataFile = list(
            "FileName" = pathParentageDataFile,
            "OutputType" = OutputType,
            "SortedBy" = OutputSortedBy,
            "IncludeNonExclusionProbabilities" = paste0(as.integer(OutputIncludeNonExclusionProbabilities))
          ))
        
        if (exists(x = "ParentageParameters", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(ParentageParameters))
        }
        if (exists(x = "OffspringFile", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(OffspringFile))
        }
        if (exists(x = "CandidateFemaleFile", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(CandidateFemaleFile))
        }
        if (exists(x = "CandidateMaleFile", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(CandidateMaleFile))
        }
        if (exists(x = "ParentageSummaryFile", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(ParentageSummaryFile))
        }
        if (exists(x = "ParentageDataFile", where = cervus_crv_file)) {
          cervus_crv_file <- within(cervus_crv_file, rm(ParentageDataFile))
        }
        
        if (!is.na(wineCommand)) {
          cervus_crv_par$OffspringFile$FileName = paste0(wineHomeDirectory, pathOffspringFile)
          cervus_crv_par$ParentageSummaryFile$FileName = paste0(wineHomeDirectory, pathParentageSummaryFile)
          cervus_crv_par$ParentageDataFile$FileName = paste0(wineHomeDirectory, pathParentageDataFile)
          
          # Normal cervus runs perfect, but currently a wine bug is breaking the command line version
          # https://bugs.winehq.org/show_bug.cgi?id=49334
          # as a work around, r will create the temp file locations needed
          # you shouldn't have to worry about temp folder bloat, as CervusCL clears it after it runs
          
          dir.create(
            path = file.path(wineTempDirectory, AnalysisFolderPath, fsep = .Platform$file.sep), 
            recursive = TRUE, 
            showWarnings = FALSE)
          
        }
        
        cervus_crv_par <- append(cervus_crv_file, cervus_crv_par)
        
        ini::write.ini(x = cervus_crv_par, filepath = pathAnalysisSettings) # requires ini package to format settings file
        
        system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/PAR /O")) # run the analysis
        system(command = paste0("cat ", '"', pathParentageSummaryFile, '"')) # display the results
        cat("\nAnalysis complete!")
        
        # No Wine run
        
        if (is.na(wineCommand)) {
          system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/PAR /O")) # run the simulation
          if (ResultsToConsole) {
            system(command = paste0("cat ", '"', pathParentageSummaryFile, '"')) # display the results
          }
        }
        
        # Wine run
        
        if (!is.na(wineCommand)) {
          system(command = paste0(wineCommand, ' "', CervusCLPath, '" ', '"', wineHomeDirectory, pathAnalysisSettings, '" ', "/PAR /O")) # run the allele frequency analysis
          if (ResultsToConsole) {
            system(command = paste0("cat ", '"', pathParentageSummaryFile, '"'))  # display the results
          }
        }
        
        cat("\nSimulations complete!\n")
        
      } else {
        cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
      }
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  }
