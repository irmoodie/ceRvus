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

SummaryStatistics$NOffspring <- read.table(text = summaryfile[(indivstatloc+2):(indivstatloc+4)], header = FALSE, sep = "\t") |> 
  tidyr::separate(col = V1, into = c("stat", "value"), sep = ":") |> 
  dplyr::mutate(stat = trimws(stat),
                value = trimws(value))

if (stringr::str_split(string = basename(PARResultsFile), pattern = "\\.")[[1]][2] == "csv") {
  SummaryStatistics$results$all <- readr::read_csv(PARResultsFile, show_col_types = FALSE)
}

if (stringr::str_split(string = basename(PARResultsFile), pattern = "\\.")[[1]][2] != "csv") {
  SummaryStatistics$results$all <- readr::read_tsv(PARResultsFile, show_col_types = FALSE)
}

SummaryStatistics$results$all <- SummaryStatistics$results$all |>
  dplyr::rename(offspring_id = "Offspring ID",
                loci_typed_offspring = "Loci typed...2",
                known_id = "Mother ID",
                loci_typed_known = "Loci typed...4",
                offspring_known_loci_compared = "Pair loci compared...5",
                offspring_known_loci_mismatching = "Pair loci mismatching...6",
                offspring_known_lod_score = "Pair LOD score...7",
                candidate_id = "Candidate father ID",
                loci_typed_candidate = "Loci typed...9",
                offspring_candidate_loci_compared = "Pair loci compared...10",
                offspring_candidate_loci_mismatching = "Pair loci mismatching...11",
                offspring_candidate_lod_score = "Pair LOD score...12",
                offspring_candidate_delta = "Pair Delta",
                offspring_candidate_confidence = "Pair confidence",
                offspring_candidate_known_loci_compared = "Trio loci compared",
                offspring_candidate_known_loci_mismatching = "Trio loci mismatching",
                offspring_candidate_known_lod_score = "Trio LOD score",
                offspring_candidate_known_delta = "Trio Delta",
                offspring_candidate_known_confidence = "Trio confidence")

SummaryStatistics$results$strict <- SummaryStatistics$results$all |> 
  dplyr::filter(offspring_candidate_known_confidence == "*")

SummaryStatistics$results$relaxed <- SummaryStatistics$results$all |> 
  dplyr::filter(offspring_candidate_known_confidence == "+" | offspring_candidate_known_confidence == "*")


return(SummaryStatistics)
}

