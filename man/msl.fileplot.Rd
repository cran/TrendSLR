% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msl.fileplot.R
\name{msl.fileplot}
\alias{msl.fileplot}
\title{Plotting to file options in JPEG format.}
\usage{
msl.fileplot(x, resol = 1800, wdir = " ", file_name = " ",
  type = 1, ci = 1, header = TRUE)
}
\arguments{
\item{x}{object of class \dQuote{msl.trend} (see \code{\link{msl.trend}}) or
\dQuote{custom.trend} (see \code{\link{custom.trend}}).}

\item{resol}{numeric, enables a user defined resolution in dpi from [300 to 1800]
where 1800 is the default setting.}

\item{wdir}{character string, providing the name of the directory to send
output files (e.g., \dQuote{C:/myproject/}). If this field is left blank the function
will terminate with a warning message sent to the console.}

\item{file_name}{is a character string indicating the name of the output file.
There is no need to include the file extension *.jpeg. If this field is left
blank the output file will be automatically saved in the defined  directory
(wdir) under the default name \dQuote{Plot1.jpeg}.}

\item{type}{numeric, enables a user defined input to select the type of chart
to be plotted. 5 seperate options are available:

\itemize{
  \item 1: The default setting provides a single 3 panel chart with the time
  series in the top panel, velocity in the middle panel and acceleration
  in the bottom panel (width = 160 mm, height = 210 mm);
  \item 2: Single panel plot of time series (width = 160 mm, height = 80 mm);
  \item 3: Single panel plot of velocity (width = 160 mm, height = 80 mm);
  \item 4: Single panel plot of acceleration (width = 160 mm, height = 80 mm); and
  \item 5: Alternative 2 panel chart with the time series in the top panel
  and velocity in the bottom panel (width = 160 mm, height = 150 mm).
}}

\item{ci}{numeric, enables a user defined input to select the type of
confidence interval to be displayed on the plots. The default setting (ci = 1)
corresponds to a 95\% confidence interval whilst ci = 2 provides a 99\%
confidence interval.}

\item{header}{logical, if \sQuote{TRUE} then the station_name (if provided)
in the \dQuote{msl.trend} (see \code{\link{msl.trend}}) or
\dQuote{custom.trend} (see \code{\link{custom.trend}}) object will be passed
to the main banner printed above the plot. If set to \sQuote{FALSE} then the
banner header will be excluded. Default = TRUE.}
}
\description{
Plotting to file options in JPEG format.
}
\details{
This function provides report quality JPEG format summary plots for both
\dQuote{msl.trend} (see \code{\link{msl.trend}}) and \dQuote{custom.trend}
(see \code{\link{custom.trend}}) objects. The same range of alternative screen
plotting options are available via \code{\link{msl.screenplot}}.
}
\examples{

# Plot to file from "custom.trend" object

data(t) # "custom.trend" object
str(t) # check object

# -------------------------------------------------------------------------
# The following call to msl.fileplot can be found in the temporary
# directory under the file name "Plot1.jpeg".
# -------------------------------------------------------------------------

wd <- tempdir() # find temp directory
msl.fileplot(t, wdir = wd) # default screen plot output

}
\seealso{
\code{\link{msl.trend}}, \code{\link{custom.trend}},
\code{\link{msl.screenplot}}
}
