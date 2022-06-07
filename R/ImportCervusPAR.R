#' Import summary and data from a Cervus parentage analysis
#'
#' Imports the summary and data files produced from a Cervus parentage analysis.
#' @import stringr
#' @import dplyr
#' @import tidyr
#' @param PARSummaryFile The summary file (.txt)
#' @param PARResultsFile The results file (.sim)
#' @return A list containing the summarised information, and the results of the parentage analysis.
#' @export
#'

ImportCervusPAR <- function(PARSummaryFile, PARResultsFile){

summaryfile <- readLines(PARSummaryFile)

sumstatloc <- grep(
  "Level       Confidence \\(%\\)  Critical Delta  Assignments        Assignment Rate",
  summaryfile
)

SummaryStatistics <- list()
for (table in sumstatloc) {
  SummaryStatistics[[paste0(stringr::str_remove(summaryfile[table-2], pattern = ":"))]] <- 
    read.table(text = summaryfile[(table+2):(table+5)],
               header = FALSE, sep = "\t")
    
}

for (table in names(SummaryStatistics)) {
  z <- data.frame()
  for (row in 1:nrow(SummaryStatistics[[table]])) {
    x <- stringr::str_remove_all(SummaryStatistics[[table]][row,], pattern = "\\(") |> 
      stringr::str_remove_all(pattern = "\\)") |> 
      stringr::str_split(pattern = "\\s+")
    
    if (length(x[[1]]) == 5) {
      y <- data.frame("Level" = x[[1]][1],
                 "Confidence" = NA,
                 "CriticalDelta" = NA,
                 "AssignmentsObserved" = x[[1]][2],
                 "AssignmentsExpected" = x[[1]][3],
                 "AssignmentRateObserved" = x[[1]][4],
                 "AssignmentRateExpected" = x[[1]][5]
                 )
    }
    
    if (length(x[[1]]) == 7) {
      y <- data.frame("Level" = x[[1]][1],
                 "Confidence" = x[[1]][2],
                 "CriticalDelta" = x[[1]][3],
                 "AssignmentsObserved" = x[[1]][4],
                 "AssignmentsExpected" = x[[1]][5],
                 "AssignmentRateObserved" = x[[1]][6],
                 "AssignmentRateExpected" = x[[1]][7]
      )
    }
    z <- rbind(z,y)
  }
  SummaryStatistics[[table]] <- z
}

indivstatloc <- grep(
  "\\*\\*\\*\\* Number of indiv",
  summaryfile
)

SummaryStatistics[["NOffspring"]] <- read.table(text = summaryfile[(indivstatloc+2):(indivstatloc+4)], header = FALSE, sep = "\t") |> 
  tidyr::separate(col = V1, into = c("stat", "value"), sep = ":") |> 
  dplyr::mutate(stat = trimws(stat),
                value = trimws(value))

SummaryStatistics[["Results"]] <- readr::read_tsv(PARResultsFile)


return(SummaryStatistics)
}

