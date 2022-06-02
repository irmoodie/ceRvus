#' Run parentage simulations using CervusCL
#'
#' You should be familiar with Cervus before using this function. 
#' The arguments relate directly to settings within Cervus, usually selected by the user using the GUI.
#' Argument names are such that it can easily be inferred what setting they change within a Cervus project file (.crv).
#' @import ini
#' @param CervusCLPath Path to cervusCL.exe (Cervus is found at http://www.fieldgenetics.com).
#' @param AnalysisFolderPath Path to the folder which contains the files required for the analysis. All output will also be saved in this folder.
#' @param AnalysisName Can specify a custom name for the analysis, which will be saved in the output filenames.
#' @param AnalysisType Set as "Maternity"  or "Paternity" or "Parent pair (sexes known)" or "Parent pair (sexes unknown)".
#' @param NOffspring Number of offspring to simulate, default to 10000, but use more for final analysis.
#' @param NCandidateMales Average number of candidate male parents per offspring, regardless of being sampled (set to 0 for maternity analysis).
#' @param PropCandidateMalesSampled Proportion of candidate male parents that are sampled (1 is not recommended, set to 0 for maternity analysis).
#' @param NCandidateFemales Average number of candidate female parents per offspring, regardless of being sampled (set to 0 for paternity analysis).
#' @param PropCandidateFemalesSampled Proportion of candidate female parents that are sampled (1 is not recommended, set to 0 for paternity analysis).
#' @param PropLociTyped Proportion of loci typed allows for missing data and should be an average value across loci and individuals. It is estimated in the ALF.
#' @param PropLociMistyped Proportion of loci mistyped allows for errors in genotyping and should be an average value across loci and individuals. See the Cervus help files for how it handles errors.
#' @param MinTypedLoci Threshold for number of typed loci for an individual to be included in the analysis.
#' @param CriticalStatisticName "Delta" or "LOD". "Delta" recommended. See the Cervus help files.
#' @param TruncateAtZero TRUE/FALSE
#' @param RelaxedConfidence Relaxed confidence level (default is 80%)
#' @param StrictConfidence Strict confidence level (default is 95%)
#' @param SimulateInbreeding TRUE/FALSE: Should inbreeding be simulated? You should also set ParentRelatedness and InbreedingRate if TRUE.
#' @param ParentRelatedness The probability that an allele in one parent is identical by descent to the corresponding allele in the other parent.
#' @param InbreedingRate The proportion of simulated offspring whose true parents are related to each other.
#' @param AlwaysTestSelfing TRUE/FALSE: Test for self-fertilisation
#' @param SimulateFemaleRelatives TRUE/FALSE: simulate relatives of the mother-father-offspring trio among the pool of candidate parents.
#' @param FemalePropRelatives The average proportion of female candidate parents that are relatives. 
#' @param FemaleRelatedTo Relatives can be related to the offspring, its true mother or its true father. See Cervus helpfiles.
#' @param FemaleRelatedness The probability that each allele belonging to a relative is identical by descent to the corresponding allele carried by the individual to which they are related. 
#' @param MalePropRelatives TRUE/FALSE: simulate relatives of the mother-father-offspring trio among the pool of candidate parents.
#' @param MaleRelatedTo Relatives can be related to the offspring, its true mother or its true father. See Cervus helpfiles.
#' @param MaleRelatedness The probability that each allele belonging to a relative is identical by descent to the corresponding allele carried by the individual to which they are related. 
#' @param UseCorrectedLikelihoods TRUE/FALSE: Use Kalinowski et al. (2007) likelihood equations (recommended)
#' @param UseMistypingRateAsLikelihoodErrorRate TRUE/FALSE
#' @param LikelihoodErrorRate Default is 0.01
#' @param RepeatSimulation TRUE/FALSE: Repeat simulation using the same parameters
#' @param NRepeats How many repeats?
#' @param ApplyPreviousSimulationData TRUE/FALSE: See Cervus helpfiles.
#' @param GenerateTables TRUE/FALSE: See Cervus helpfiles.
#' @param SaveRawStatisticScores TRUE/FALSE: See Cervus helpfiles.
#' @param GenerateHistograms TRUE/FALSE: See Cervus helpfiles.
#' @param NCategories Histogram setup
#' @param MinStatistic Histogram setup
#' @param MaxStatistic Histogram setup
#' @param PreviousSimulationDataFile Path to previous simulation data file (not required unless ApplyPreviousSimulationData = TRUE)
#' @return Currently prints the results to the console, and saves the simulation data in AnalysisFolderPath
#'
#' @export

CervusSIM <- 
  function(
    CervusCLPath,
    AnalysisFolderPath, 
    AnalysisName = "Cervus_Analysis",
    AnalysisType,
    NOffspring = 10000,
    NCandidateMales = 0,
    PropCandidateMalesSampled = 0,
    NCandidateFemales = 0,
    PropCandidateFemalesSampled = 0,
    PropLociTyped,
    PropLociMistyped = 0.01,
    MinTypedLoci,
    CriticalStatisticName = "Delta",
    TruncateAtZero = TRUE,
    RelaxedConfidence = 80,
    StrictConfidence = 95,
    SimulateInbreeding = FALSE,
    ParentRelatedness = 0,
    InbreedingRate = 0,
    AlwaysTestSelfing = FALSE,
    SimulateFemaleRelatives = FALSE,
    FemalePropRelatives = 0,
    FemaleRelatedTo = "Offspring",
    FemaleRelatedness = 0,
    SimulateMaleRelatives = FALSE,
    MalePropRelatives = 0,
    MaleRelatedTo = "Offspring",
    MaleRelatedness = 0,
    UseCorrectedLikelihoods = TRUE,
    UseMistypingRateAsLikelihoodErrorRate = TRUE,
    LikelihoodErrorRate = 0.01,
    RepeatSimulation = FALSE,
    NRepeats = 1,
    ApplyPreviousSimulationData = FALSE,
    GenerateTables = FALSE,
    SaveRawStatisticScores = FALSE,
    GenerateHistograms = FALSE,
    NCategories = 5,
    MinStatistic = 0,
    MaxStatistic = 0,
    PreviousSimulationDataFile = ""){
    
    if (!missing(CervusCLPath)) {
      if (file.exists(CervusCLPath)) {
        
        pathAnalysisSettings <- file.path(AnalysisFolderPath, paste0(AnalysisName,"_settings", ".crv"), fsep = "\\")
        pathSimulationSummaryFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Simulation.txt"), fsep = "\\")
        pathSimulationDataFile <- file.path(AnalysisFolderPath, paste0(AnalysisName, "_Simulation.sim"), fsep = "\\")
        
        CervusCRVFile <- ini::read.ini(pathAnalysisSettings)
        
        CervusCRVSIM <- list(
          SimulationParameters = list(
            "AnalysisType" = paste0(AnalysisType),
            "NOffspring" = paste0(NOffspring),
            "NCandidateMales" = paste0(NCandidateMales),
            "PropCandidateMalesSampled" = paste0(PropCandidateMalesSampled),
            "NCandidateFemales" = paste0(NCandidateFemales),
            "PropCandidateFemalesSampled" = paste0(PropCandidateFemalesSampled),
            "PropLociTyped" = paste0(PropLociTyped),
            "PropLociMistyped" = paste0(PropLociMistyped),
            "MinTypedLoci" = paste0(MinTypedLoci),
            "CriticalStatisticName" = paste0(CriticalStatisticName),
            "TruncateAtZero" = paste0(as.integer(TruncateAtZero)),
            "RelaxedConfidence" = paste0(RelaxedConfidence),
            "StrictConfidence" = paste0(StrictConfidence),
            "SimulateInbreeding" = paste0(as.integer(SimulateInbreeding)),
            "ParentRelatedness" = paste0(ParentRelatedness),
            "InbreedingRate" = paste0(InbreedingRate),
            "AlwaysTestSelfing" = paste0(as.integer(AlwaysTestSelfing)),
            "SimulateFemaleRelatives" = paste0(as.integer(SimulateFemaleRelatives)),
            "FemalePropRelatives" = paste0(FemalePropRelatives),
            "FemaleRelatedTo" = paste0(FemaleRelatedTo),
            "FemaleRelatedness" = paste0(FemaleRelatedness),
            "SimulateMaleRelatives" = paste0(as.integer(SimulateMaleRelatives)),
            "MalePropRelatives" = paste0(MalePropRelatives),
            "MaleRelatedTo" = paste0(MaleRelatedTo),
            "MaleRelatedness" = paste0(MaleRelatedness),
            "UseCorrectedLikelihoods" = paste0(as.integer(UseCorrectedLikelihoods)),
            "UseMistypingRateAsLikelihoodErrorRate" = paste0(as.integer(UseMistypingRateAsLikelihoodErrorRate)),
            "LikelihoodErrorRate" = paste0(LikelihoodErrorRate)
          ),
          SimulationSummaryFile = list(
            "FileName" = pathSimulationSummaryFile # must be absolute path (e.g. C:\\...)
          ),
          SimulationDataFile = list(
            "FileName" = pathSimulationDataFile # must be absolute path (e.g. C:\\...)
          ),
          SimulationOutput = list(
            "RepeatSimulation" = paste0(as.integer(RepeatSimulation)),
            "NRepeats" = paste0(NRepeats),
            "ApplyPreviousSimulationData" = paste0(as.integer(ApplyPreviousSimulationData)),
            "GenerateTables" = paste0(as.integer(GenerateTables)),
            "SaveRawStatisticScores" = paste0(as.integer(SaveRawStatisticScores)),
            "GenerateHistograms" = paste0(as.integer(GenerateHistograms)),
            "NCategories" = paste0(NCategories),
            "MinStatistic" = paste0(MinStatistic),
            "MaxStatistic" = paste0(MaxStatistic)
          ),
          PreviousSimulationDataFile = list(
            "FileName" = paste0(PreviousSimulationDataFile) # must be absolute path (e.g. C:\\...)
          )
        )
        
        if (exists(x = "SimulationParameters", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(SimulationParameters))
        }
        if (exists(x = "SimulationSummaryFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(SimulationSummaryFile))
        }
        if (exists(x = "SimulationDataFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(SimulationDataFile))
        }
        if (exists(x = "SimulationOutput", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(SimulationOutput))
        }
        if (exists(x = "PreviousSimulationDataFile", where = CervusCRVFile)) {
          CervusCRVFile <- within(CervusCRVFile, rm(PreviousSimulationDataFile))
        }
        
        CervusCRVSIM <- append(CervusCRVFile, CervusCRVSIM)
        
        ini::write.ini(x = CervusCRVSIM, filepath = pathAnalysisSettings) # requires ini package to format settings file
        
        system(command = paste0('"', CervusCLPath, '" ', '"', pathAnalysisSettings, '" ', "/SIM /O")) # run the simulation
        system(command = paste0("cat ", '"', pathSimulationSummaryFile, '"')) # display the results
        cat("\nSimulations complete!")
        
      } else {
        cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
      }
    } else {
      cat("WARNING: Cannot locate CervusCL.exe.\nPlease ensure you have provided the full system path e.g. C:\\...")
    }
    
  }
