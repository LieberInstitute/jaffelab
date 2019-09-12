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
#' @param sheet_name A `character(1)` vector with the name of the Google
#' Drive spreadsheet that will be updated (or created if it doesn't exist).
#' @param n The number of presenters per meeting.
#' @param repeat_day An `integer(1)` with the number of days between lab
#' meetings. Use 7 for a weekly lab meeting.
#'
#' @return A one row table with the `googledrive` information for the file
#' that was uploaded.
#' @export
#' @importFrom googledrive drive_find drive_upload drive_update
#'
#' @examples
#'
#' ## Presenters as of 2019-09-12
#' jaffe_research_presenters <- c('Brianna', 'Emily', 'Josh', 'Kira', 'Leo',
#'     'Maddy', 'Matt', 'Nick')
#'
#' \dontrun{
#' ## You'll need to have access to Google Drive through
#' ## googledrive::drive_auth() set up.
#'
#'
#' ## Update lab presenters sheet
#' lab_presenters(
#'     presenters = jaffe_research_presenters,
#'     start_date = '2019-09-19',
#'     sheet_name = 'Jaffelab research presenters',
#'     n = 2
#' )
#' }
#'

lab_presenters <- function(presenters, start_date = '2019-09-19',
    sheet_name = 'Jaffelab research presenters', n = 2, repeat_day = 7) {
    ## Check inputs
    if(!is.character(presenters))
        stop("'presenters' should be a character vector.", call. = FALSE)
    if(!is.character(start_date))
        stop("'start_date' should be a character vector.", call. = FALSE)
    if(!is.character(sheet_name))
        stop("'sheet_name' should be a character vector.", call. = FALSE)
    if(n < 1) stop("'n' should be at least 1.", call. = FALSE)
    if(!identical(nchar(strsplit(start_date, '-')[[1]]), c(4L, 2L, 2L)))
        stop("'start_date' should be in the format YYYY-MM-DD.", call. = FALSE)

    ## Keep the unique and sort them to make it 100% reproducible
    presenters <- sort(unique(presenters))
    date_start <- as.Date(start_date)

    ## Define groups of presenters
    i <- rep(seq(from = 1, to = length(presenters), by = n), each = n)

    ## Randomize presenters
    set.seed(as.numeric(gsub('-', '', start_date)))
    presenters <- sample(presenters)

    ## Group presenters (add NAs if missing)
    if(length(i) > length(presenters)) {
        presenters <- c(presenters, rep(NA, length(i) - length(presenters)))
    }
    presenters_grouped <- split(presenters, i)

    df <- data.frame(do.call(rbind, presenters_grouped), stringsAsFactors = FALSE)
    colnames(df) <- paste0('Presenter ', seq_len(n))
    df$Date <- date_start + repeat_day * (seq_len(length(presenters_grouped)) - 1)

    ## Write temp csv
    tmp_csv <- file.path(tempdir(), paste0(sheet_name, '.csv'))
    write.csv(df, file = tmp_csv, row.names = FALSE)

    ## Find sheets with the name
    gdoc <- googledrive::drive_find(
        pattern = paste0('^', sheet_name, '$'),
        type = "spreadsheet"
    )

    if(nrow(gdoc) == 0) {
        result <- googledrive::drive_upload(tmp_csv, path = ,
            type = 'spreadsheet', name = sheet_name)
    } else if (nrow(gdoc) > 1) {
        stop("The 'sheet_name' you chose is not unique.\n",
            "Check your Google Drive files!", .call = FALSE)
    } else {
        result <- googledrive::drive_update(gdoc, tmp_csv)
    }

    return(result)
}
