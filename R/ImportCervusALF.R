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
    header = TRUE, na.strings = "ND")
  
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

  ReadPerLocusAlleleFrequencies <- function(x, k, r, col_names){
    text <- x[r[-(1:3)][1:k]]
    text <- gsub("Not done", "ND", text)
    return(read.table(text = text, col.names = col_names, na.strings = "ND"))
  }
  
  for (loci_n in 1:length(x)) {
    
    LociName <- LociSummary$Locus[loci_n]
    k_alleles <- LociSummary$k[loci_n]
    cols <- c("Allele", "Count", "Hets", "Homs", "Freq", "FreqWithNull")

    # extract the allele frequency table for each loci
    if (loci_n != length(x)) {
      range <- (x[[loci_n]]):(x[[loci_n+1]])
    }
    if (loci_n == length(x)) {
      range <- (x[[loci_n]]):(grep("\\*\\*\\*\\*\\*\\*\\*\\*", summaryfile))
    }
    if (k_alleles > 0) {
      PerLocusAlleleFrequencies[[LociName]] <- 
        ReadPerLocusAlleleFrequencies(x = summaryfile, r = range, k = k_alleles, col_names = cols)
    }
    if (k_alleles == 0) {
      PerLocusAlleleFrequencies[[LociName]] <- data.frame(matrix(NA, ncol = length(cols), nrow = 1))
      colnames(PerLocusAlleleFrequencies[[LociName]]) <- cols
    }
    
    # extract the statistics table for each loci
    stats_start <- 3 + k_alleles + 2
    null_allele_line <- grep("Null allele frequency estimate", summaryfile[range])
    stats_lines <- summaryfile[range[stats_start:null_allele_line]]
    
    PerLocusStatistics[[LociName]] <- stringr::str_split(
      string = stats_lines,
      pattern = ":",
      simplify = TRUE)
    PerLocusStatistics[[LociName]][,2] <- stringr::str_remove_all(string = PerLocusStatistics[[LociName]][,2], pattern = " ")
    PerLocusStatistics[[LociName]][,2] <- stringr::str_replace_all(string = PerLocusStatistics[[LociName]][,2], pattern = "Notdone", replacement = NA_character_)
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
