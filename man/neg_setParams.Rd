% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neg_setParams.R
\name{neg_setParams}
\alias{neg_setParams}
\title{Optimize parameters for predicting negative voice onset time}
\usage{
neg_setParams(
  directory = ".",
  when_use_f0 = 10,
  verbose = TRUE,
  plot = TRUE,
  ...
)
}
\arguments{
\item{directory}{A string giving the location of a directory with WAV/TextGrid
pairs. Ideally, each WAV file should contain just one syllable with a stop,
where the stop closure is placed somewhere in the first half of the sound file.
The TextGrid tier containing hand-segmented VOT should be named \code{vot}.
Default is the current working directory.}

\item{when_use_f0}{Numeric; default is \code{10}. Per default, voicing onset is
predicted with \link{negativeVOT} using \code{method='acf'}. This is the preferred method
because it is significantly faster than the alternative \code{method='f0'}. For
this reason, parameters with \code{method='f0'} are only tested if the average
prediction error with \code{method='acf'} is above a certain threshold.
\code{when_use_f0} sets this threshold; i.e. by default, parameters for
\code{method='f0'} are tested when the average prediction error with \code{method='acf'}
is above 10 ms.}

\item{verbose}{Logical; should messages be printed in the console indicating
the progress? Default is \code{TRUE}.}

\item{plot}{Logical; should predicted VOT using the resulting parameters
be plotted against the hand-annotated VOT?}

\item{...}{Additional graphical parameters to be passed on to
\link{plot_training_diff}.}
}
\value{
A named list with optimized parameters for the training data; see
\link{negativeVOT} for more information.
}
\description{
If the user provides a few examples of WAV/TextGrid pairs with hand-annotated
voice onset time of stops, this function will attempt to predict VOT with
a large range of parameters, and select those parameters that on average
minimize prediction error the most. These parameters can then be used for
the rest of the data set.
}
\examples{
datapath <- system.file('extdata/vd', package='getVOT')
neg_setParams(directory=datapath)
}
\seealso{
The returned list can be passed to the \link{negativeVOT} argument
\code{params_list}, or the argument \code{neg_params_list} in other functions
\link{getVOT}, \link{VOT2newTG}, \link{addVOT2TG}, \link{addVOT2emuDB}. \link{pos_setParams} is a
sister function for optimizing parameters for predicting positive voice
onset time. \link{plot_training_diff} can be used to illustrate how closely the
optimized parameters predict the hand-annotated training data.
}
