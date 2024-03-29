% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lab_presenters.R
\name{lab_presenters}
\alias{lab_presenters}
\title{Randomize presenters and upload them to Google Drive}
\usage{
lab_presenters(
  presenters,
  start_date = "2019-09-18",
  sheet_name = "Jaffelab research presenters",
  n = 2,
  repeat_day = 7
)
}
\arguments{
\item{presenters}{A \code{character()} vector of presenter names.}

\item{start_date}{A \code{character(1)} vector with the start date in a
\code{YYYY-MM-DD} format.}

\item{sheet_name}{A \code{character(1)} vector with the name of the Google
Drive spreadsheet that will be updated (or created if it doesn't exist).}

\item{n}{The number of presenters per meeting.}

\item{repeat_day}{An \code{integer(1)} with the number of days between lab
meetings. Use 7 for a weekly lab meeting.}
}
\value{
A one row table with the \code{googledrive} information for the file
that was uploaded.
}
\description{
This function takes a set of presenter names and randomizes them into groups
of size \code{n} based on the \code{start_date}. The resulting table is then uploaded
to Google Drive. If the \code{sheet_name} exists, it's updated, otherwise a new
spreadsheet is created which you can share later.
}
\examples{

## Presenters as of 2019-09-12
jaffe_research_presenters <- c(
    "Brianna", "Emily", "Josh", "Kira", "Leo",
    "Maddy", "Matt", "Nick"
)

if (googledrive::drive_has_token()) {
    ## You'll need to have access to Google Drive through
    ## googledrive::drive_auth() set up.

    ## Set the seed for reproducibility of the results
    set.seed(20190918)

    ## Update lab presenters sheet
    lab_presenters(
        presenters = jaffe_research_presenters,
        start_date = "2019-09-18",
        sheet_name = "Jaffelab research presenters - example",
        n = 2
    )

    ## If it's a new sheet, we recommend sharing an editable link so other
    ## lab members can swap out as necessary.
}
}
\seealso{
googledrive_csv
}
\author{
Leonardo Collado-Torres
}
