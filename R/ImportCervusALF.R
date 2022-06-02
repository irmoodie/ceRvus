#' Import summary of allele frequency analysis
#'
#' Imports the summary file produced from a Cervus allele frequency analysis (the .txt file). Can be called directly in the CervusALF function with ImportALF = TRUE.
#' @import stringr
#' @param CervusALFSummaryFile The path to the allele frequency analysis summary file (.txt file)
#' @return Returns the tables in the Cervus allele frequency analysis summary file, formatted for use in R
#' @examples 
#' ImportCervusALF(ALFSummaryFile = "AlleleFrequencyAnalysis.txt")
#' @export

ImportCervusALF <- function(ALFSummaryFile){
  
  summaryfile <- readLines(ALFSummaryFile)
  
  # Isolate the summary statistics and coerce into a dataframe
  LociSummary <- read.table(
    text = summaryfile[(grep("\\*\\*\\*\\* Summary statistics \\*\\*\\*\\*", summaryfile)+3):(grep("Number of individuals:", summaryfile)-2)],
    header = TRUE)
  
  # Isolate the dataset statistics and coerce into a dataframe
  DataSummary <- stringr::str_split(
    string = summaryfile[(grep("Number of individuals:", summaryfile)):grep("Combined non-exclusion probability \\(sib identity\\):", summaryfile)],
    pattern = ":",
    simplify = TRUE)
  DataSummary[,2] <- stringr::str_remove_all(string = DataSummary[,2], pattern = " ")
  DataSummary <- data.frame(stat = DataSummary[,1], value = DataSummary[,2])
  
  # Find the line where each loci stats start
  x <- vector()
  for (loci in LociSummary$Locus) {
    x[loci] <- grep(paste0("\\*\\*\\*\\* Locus ", loci, " \\*\\*\\*\\*"), summaryfile)
  }
  
  # create lists for per locus stats
  PerLocusAlleleFrequencies <- list()
  PerLocusStatistics <- list()
  
  for (loci_n in 1:length(x)) {
    
    LociName <- LociSummary$Locus[loci_n]
    cols <- c("Allele", "Count", "Hets", "Homs", "Freq", "FreqWithNull")
    
    # extract the allele frequency table for each loci
    if (loci_n != length(x)) {
      range <- (x[[loci_n]]):(x[[loci_n+1]])
      PerLocusAlleleFrequencies[[LociName]] <- read.table(text = summaryfile[head(range[-(1:4)], -23)], col.names = cols)
    }
    
    # extract the allele frequency table for the final loci in the dataset
    if (loci_n == length(x)) {
      range <- (x[[loci_n]]):(grep("\\*\\*\\*\\*\\*\\*\\*\\*", summaryfile))
      PerLocusAlleleFrequencies[[LociName]] <- read.table(text = summaryfile[head(range[-(1:4)], -23)], col.names = cols)
    }
    
    # extract the statistics table for each loci
    PerLocusStatistics[[LociName]] <- stringr::str_split(
      string = summaryfile[head(range[-(1:(length(head(range[-(1:4)], -23))+5))], -3)],
      pattern = ":",
      simplify = TRUE)
    PerLocusStatistics[[LociName]][,2] <- stringr::str_remove_all(string = PerLocusStatistics[[LociName]][,2], pattern = " ")
    PerLocusStatistics[[LociName]] <- data.frame(stat = PerLocusStatistics[[LociName]][,1], value = PerLocusStatistics[[LociName]][,2])
    
  }
  
  # Summarise into a single object
  CervusALFSummary <- list(
    "SummaryStatistics" = list(
      "LociSummary" = LociSummary,
      "DataSummary" = DataSummary
    ),
    "PerLocusStatistics" = list(
      "PerLocusAlleleFrequencies" = PerLocusAlleleFrequencies,
      "PerLocusStatistics" = PerLocusStatistics
    )
  )
  
  return(CervusALFSummary)
}
