#' Upload a CSV file to googledrive
#'
#' This function uploads a data.frame as a CSV file to Google Drive. If the
#' name of the spreadsheet exists already, then it gets updated. Otherwise
#' it uploads it for the first time. Unlike using
#' `drive_upload(overwrite = TRUE)` directly, this function allows you to
#' keep the sharing and publishing settings of a spreadsheet if it exists
#' already.
#'
#' @param df A `data.frame()` to upload.
#' @param sheet_name A `character(1)` vector with the name of the Google
#' Drive spreadsheet that will be updated (or created if it doesn't exist).
#'
#' @return  A one row table with the `googledrive` information for the file
#' that was uploaded.
#' @export
#' #' @importFrom googledrive drive_find drive_upload drive_update
#' @author Leonardo Collado-Torres
#'
#' @examples
#'
#' if(googledrive::drive_has_token()) {
#' ## You'll need to have access to Google Drive through
#' ## googledrive::drive_auth() set up.
#'
#'
#' ## Upload a table to google drive
#' googledrive_csv(mtcars, paste(Sys.Date(), 'jaffelab::googledrive_csv() example'))
#'
#' }
#'
googledrive_csv <- function(df, sheet_name) {

    if(!is.character(sheet_name))
        stop("'sheet_name' should be a character vector.", call. = FALSE)

    ## Write temp csv
    tmp_csv <- file.path(tempdir(), paste0(sheet_name, '.csv'))
    write.csv(df, file = tmp_csv, row.names = FALSE)

    ## Find sheets with the name
    gdoc <- googledrive::drive_find(
        pattern = paste0('^', sheet_name, '$'),
        type = "spreadsheet"
    )

    ## Upload if absent, update if present
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
