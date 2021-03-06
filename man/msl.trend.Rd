% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msltrend.R
\name{msl.trend}
\alias{msl.trend}
\title{Isolate trend component from mean sea level records.}
\usage{
msl.trend(object, station_name = " ", fillgaps = 1, iter = 10000,
  vlm = " ", plot = TRUE, wdir = " ", save_summary = "TRUE")
}
\arguments{
\item{object}{an annual average mean sea level time series (refer \code{\link[stats]{ts}})
with water levels (in millimetres). Missing data and maximum missing data gap
are limited to 15\% and 5\%, respectively, of the data record.\cr


\strong{Warning: } If input data files do not conform to these pre-conditions,
the analysis will be terminated. It should be further noted that the existence
of long period oscillations in global mean sea level have been well recognised
in the literature (eg. Chambers et al. (2012); Minobe (1999)). Therefore, in
order to be effective for climate change and sea level research, time series
input files are recommended to have a minimum length of at least 80 years
in order that the package can identify and isloate such signals. Time series
less than 80 years in length will be analysed but a warning will be displayed.}

\item{station_name}{character string, providing the name of the data record.

\strong{Note: }This field can be left blank, however, it is retained for use
in banner labelling of plotting and associated outputs.}

\item{fillgaps}{numeric, provides 5 alternative gap filling procedures for
missing data. The following options are available:

\itemize{
  \item 1: The default procedure is based on iterative gap filling using Singular Spectrum Analysis (refer \code{\link[Rssa]{igapfill}});
  \item 2: linear interpolation (refer \code{\link[zoo]{na.approx}});
  \item 3: Cubic spline interpolation (refer \code{\link[zoo]{na.approx}});
  \item 4: Stineman's interpolation (refer \code{\link[imputeTS]{na.interpolation}}); and
  \item 5: Weighted moving average (refer \code{\link[imputeTS]{na.ma}}).
}

\strong{Note: }Gap filled portions of the time series are denoted in red on
the default screen plot. This is done specifically to provide ready visual
observation to discern if the selected gap filling method provides an
appropriate estimate within the gaps in keeping with the remainder of the
historical record. Depending on the nature of the record and extent of gaps,
some trial and error between alternatives might be necessary to optimise gap
filling. This is best achieved via the \code{\link{gap.fillview}} function
which permits visual screen checking of various gap-filling options prior to
undertaking the full trend analysis.}

\item{iter}{numeric, enables a user defined number of iterations for
bootstrapping to determine error margins. The user range is [500 to 10000]
where 10000 is the default setting.\cr

\strong{Warning: }Although the default setting provides a more accurate basis
for estimating error margins, the degree of iterations slows the analysis and
can take several minutes to run.}

\item{vlm}{numeric, enables a user defined quantum for vertical land motion
in mm/year within the range [-20 to 20]. This rate is used to convert the rate
of relative sea level rise to an estimate of geocentric sea level rise. Positive
rates of vlm are associated with land uplift, while conversely negative rates
of vlm are associated with subsidence. This can be left blank in which case
only estimates of relative mean sea level will be determined.}

\item{plot}{logical, if \dQuote{TRUE} then the original time series is
plotted to the screen along with the trend component and the result of gap
filling (where necessary). 95\% confidence intervals have also been applied.
Default = \dQuote{TRUE}.}

\item{wdir}{character string, providing the name of the directory to send
output files (e.g., \dQuote{C:/myproject/}) when the save_summary argument is set to "TRUE".
If this field is left blank the save_summary argument is switched off and a message
will be sent to the console.}

\item{save_summary}{logical, if \dQuote{TRUE} the object$Summary portion
of the returned value is exported direct to the defined  directory (wdir) and
saved as "detailed_summary_output.csv". Default = \dQuote{FALSE}.}
}
\value{
An object of class \dQuote{msl.trend} is returned with the following
elements:
\describe{
 \item{\strong{$Station.Name: }}{the name of the data record.}
 \item{\strong{$Summary: }}{a summary data frame of the relevant attributes
 relating to the trend and
  the inputted annual average data set, including:}
 \itemize{
  \item{$Year: input data; }
  \item{$MSL: input data; }
  \item{$Trend: mean sea level trend; }
  \item{$TrendSD: standard deviation of the determined mean sea level
  trend; }
  \item{$Vel: relative velocity (or first derivative) of mean sea level trend
  (mm/year); }
  \item{$VelSD: standard deviation of the velocity of the mean sea level
  trend; }
  \item{$Acc: acceleration (or second derivative) of mean sea level trend
  (mm/year/year); }
  \item{$AccSD: standard deviation of the acceleration of the mean sea level
  trend; }
  \item{$Resids: time series of uncorrelated residuals; and }
  \item{$FilledTS: gap-filled time series (where necessary). }
  \item{$VelGeo: geocentric velocity (or first derivative) of mean sea level
  trend (mm/year)(only where vertical land motion has been supplied). }
   }
}

\describe{
 \item{\strong{$Relative.Velocity: }}{outputs the peak relative velocity
 and the year in which it occurred.}
 \item{\strong{$Vertical.Land.Motion: }}{outputs the vertical land motion
 used to convert relative to geocentric velocity (user supplied input).}
 \item{\strong{$Geocentric.Velocity: }}{outputs the peak geocentric velocity
 and the year in which it occurred (if vertical land motion supplied).}
 \item{\strong{$Acceleration: }}{outputs the peak acceleration and the
 year in which it occurred.}
 \item{\strong{$Record.Length: }}{outputs details of the start, end and
 length of the input data set.}
 \item{\strong{$Fillgaps: }}{outputs the extent of missing data (years) in
 the original record and the gap filling method used (where necessary).}
 \item{\strong{$Bootstrapping.Iterations: }}{outputs the number of iterations
 used to generate the respective standard deviations for error margins.}
 \item{\strong{$Changepoints: }}{outputs the number and time at which
 changepoints in the variance of the uncorrelated residuals occur (if any).
 Where changepoints are identified, block bootstrapping procedures are used
 with residuals quarantined between changepoints.}
   }
}
\description{
Isolate trend component from mean sea level records.
}
\details{
This function deconstructs annual average time series data into a trend
and associated velocities and accelerations, filling necessary internal structures
to facilitate all other functions in this package. The trend is isloated using
Singular Spectrum Analysis, in particular, aggregating components whose low
frequency band [0 to 0.01] exceed a threshold contribution of 75\%. Associated
velocities and accelerations are determined through the fitting of a cubic
smoothing spline to the trend with 1 degree of freedom per every 8 years of
record length. The fixed settings built into this function are based on the
detailed research and development summarised in Watson (2016a,b; 2018).
}
\examples{

data(Balt) # Baltimore mean sea level record
ts1 <- ts(Balt[2], start = Balt[1, 1]) # create time series input object
\donttest{s <- msl.trend(ts1, fillgaps = 3, iter = 500, 'BALTIMORE, USA')}

data(s)
str(s) # check structure of object
msl.screenplot(s) # check screen output

}
\references{
Chambers, D.P., Merrifield, M.A., and Nerem, R.S., 2012. Is there
a 60 year oscillation in global mean sea level? \emph{Geophysical Research Letters}, 39(18).

Minobe, S., 1999. Resonance in bidecadal and pentadecadal climate oscillations
over the North Pacific: Role in climatic regime shifts. \emph{Geophysical Research Letters},
26(7), pp.855-858.

Watson, P.J., 2016a. Identifying the best performing time series
analytics for sea-level research. In: \emph{Time Series Analysis and
Forecasting, Contributions to Statistics}, pp. 261-278, ISBN 978-3-319-28725-6.
Springer International Publishing.

Watson, P.J., 2016b. How to improve estimates of real-time acceleration in
the mean sea level signal. In: Vila-Concejo, A., Bruce, E., Kennedy, D.M.,
and McCarroll, R.J. (eds.), Proceedings of the 14th International Coastal
Symposium (Sydney, Australia). \emph{Journal of Coastal Research},
Special Issue, No. 75, pp. 780-785. Coconut Creek (Florida), ISSN 0749-0208.

Watson, P.J., 2018. \emph{Improved Techniques to Estimate Mean Sea Level,
Velocity and Acceleration from Long Ocean Water Level Time Series to Augment
Sea Level (and Climate Change) Research.} PhD Thesis, University of New South
Wales, Sydney, Australia.
}
\seealso{
\code{\link{custom.trend}}, \code{\link{gap.fillview}}, \code{\link{check.decomp}},
\code{\link{s}}, \code{\link[stats]{ts}}, \code{\link{msl.fileplot}},
\code{\link{msl.screenplot}}, \code{\link{summary}}, \code{\link{Balt}},
\code{\link[zoo]{na.approx}}, \code{\link[imputeTS]{na.interpolation}},
\code{\link[imputeTS]{na.ma}}.
}
