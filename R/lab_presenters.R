#' Randomize presenters and upload them to Google Drive
#'
#' This function takes a set of presenter names and randomizes them into groups
#' of size `n` based on the `start_date`. The resulting table is then uploaded
#' to Google Drive. If the `sheet_name` exists, it's updated, otherwise a new
#' spreadsheet is created which you can share later.
#'
#' @param presenters A `character()` vector of presenter names.
#' @param start_date A `character(1)` vector with the start date in a
#' `YYYY-MM-DD` format.
#' @param n The number of presenters per meeting.
#' @param repeat_day An `integer(1)` with the number of days between lab
#' meetings. Use 7 for a weekly lab meeting.
#' @inheritParams googledrive_csv
#'
#' @return A one row table with the `googledrive` information for the file
#' that was uploaded.
#' @export
#' @author Leonardo Collado-Torres
#' @seealso googledrive_csv
#'
#' @examples
#'
#' ## Presenters as of 2019-09-12
#' jaffe_research_presenters <- c(
#'     "Brianna", "Emily", "Josh", "Kira", "Leo",
#'     "Maddy", "Matt", "Nick"
#' )
#'
#' if (googledrive::drive_has_token()) {
#'     ## You'll need to have access to Google Drive through
#'     ## googledrive::drive_auth() set up.
#'
#'     ## Set the seed for reproducibility of the results
#'     set.seed(20190918)
#'
#'     ## Update lab presenters sheet
#'     lab_presenters(
#'         presenters = jaffe_research_presenters,
#'         start_date = "2019-09-18",
#'         sheet_name = "Jaffelab research presenters - example",
#'         n = 2
#'     )
#'
#'     ## If it's a new sheet, we recommend sharing an editable link so other
#'     ## lab members can swap out as necessary.
#' }
lab_presenters <- function(presenters, start_date = "2019-09-18",
    sheet_name = "Jaffelab research presenters", n = 2, repeat_day = 7) {
    ## Check inputs
    if (!is.character(presenters)) {
        stop("'presenters' should be a character vector.", call. = FALSE)
    }
    if (!is.character(start_date)) {
        stop("'start_date' should be a character vector.", call. = FALSE)
    }
    if (n < 1) stop("'n' should be at least 1.", call. = FALSE)
    if (!identical(nchar(strsplit(start_date, "-")[[1]]), c(4L, 2L, 2L))) {
        stop("'start_date' should be in the format YYYY-MM-DD.", call. = FALSE)
    }

    ## Keep the unique and sort them to make it 100% reproducible
    presenters <- sort(unique(presenters))
    date_start <- as.Date(start_date)

    ## Define groups of presenters
    i <- rep(seq(from = 1, to = length(presenters), by = n), each = n)

    ## Randomize presenters
    presenters <- sample(presenters)

    ## Group presenters (add NAs if missing)
    if (length(i) > length(presenters)) {
        presenters <- c(presenters, rep(NA, length(i) - length(presenters)))
    }
    presenters_grouped <- split(presenters, i)

    df <- data.frame(do.call(rbind, presenters_grouped), stringsAsFactors = FALSE)
    colnames(df) <- paste0("Presenter ", seq_len(n))
    df$Date <- date_start + repeat_day * (seq_len(length(presenters_grouped)) - 1)

    googledrive_csv(df = df, sheet_name = sheet_name)
}
