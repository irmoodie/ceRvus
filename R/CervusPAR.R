#' Run parentage analysis using CervusCL
#'
#' You should be familiar with Cervus before using this function.
#' The arguments relate directly to settings within Cervus, usually selected by the user using the GUI.
#' Argument names are such that it can easily be inferred what setting they change within a Cervus project file (.crv).
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
#' @param OutputType Default is "The most-likely parent", see the Cervus help files for more. For fractional assignment, use "All parents".
#' @param OutputSortedBy Default is "ID", see the Cervus help files for more.
#' @param OutputIncludeNonExclusionProbabilities TRUE/FALSE
#' @return Currently prints the results to the console, and saves the parentage analysis data in AnalysisFolderPath
#'
#' @export

CervusPAR <-
  function(CervusCLPath,
           AnalysisFolderPath,
           AnalysisName = "Cervus_Analysis",
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
           OutputSortedBy = "ID",
           OutputIncludeNonExclusionProbabilities = FALSE) {
    if (!missing(CervusCLPath)) {
      if (file.exists(CervusCLPath)) {
        pathAnalysisSettings <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_settings", ".crv"), fsep = "\\")
        pathParentageSummaryFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Parentage.txt"), fsep = "\\")
        pathParentageDataFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Parentage.sim"), fsep = "\\")

        pathOffspringFile <- file.path(AnalysisFolderPath, OffspringFile_FileName, fsep = "\\")

        if (CandidateFemaleFile_FileName != "") {
          pathCandidateFemaleFile <- file.path(AnalysisFolderPath, CandidateFemaleFile_FileName, fsep = "\\")
        } else {
          pathCandidateFemaleFile <- CandidateFemaleFile_FileName
        }

        if (CandidateMaleFile_FileName != "") {
          pathCandidateMaleFile <- file.path(AnalysisFolderPath, CandidateMaleFile_FileName, fsep = "\\")
        } else {
          pathCandidateMaleFile <- CandidateMaleFile_FileName
        }

        CervusCRVFile <- ini::read.ini(pathAnalysisSettings)

        CervusCRVPAR <- list(
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
          )
        )

        if (exists(x = "ParentageParameters", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(ParentageParameters))
        }
        if (exists(x = "OffspringFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(OffspringFile))
        }
        if (exists(x = "CandidateFemaleFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(CandidateFemaleFile))
        }
        if (exists(x = "CandidateMaleFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(CandidateMaleFile))
        }
        if (exists(x = "ParentageSummaryFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(ParentageSummaryFile))
        }
        if (exists(x = "ParentageDataFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(ParentageDataFile))
        }

        CervusCRVPAR <- append(CervusCRVFile, CervusCRVPAR)

        ini::write.ini(x = CervusCRVPAR, filepath = pathAnalysisSettings) # requires ini package to format settings file

        system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/PAR /O")) # run the analysis
        system(command = paste0("cat ", '"', pathParentageSummaryFile, '"')) # display the results
        cat("\nAnalysis complete!")
      } else {
        cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
      }
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
  }
