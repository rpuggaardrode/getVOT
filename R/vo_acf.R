#' Predict voicing onset in a stop based on autocorrelation
#'
#' Predicts the location of voicing onset in a short sound file from the mean
#' degree of energy autocorrelation in short intervals of the signal relative to
#' the mostautocorrelated interval during the vowel.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param start Numeric value giving the first sample of the search space in
#' the sound signal.
#' @param end Numeric value giving the final sample of the search space in
#' the sound signal.
#' @param baseline Numeric value giving the sample of a preceding landmark in
#' the sound signal (predicted closure when predicting negative VOT,
#' predicted release when predicting positive VOT).
#' @param step The duration of each chunk in terms of samples (corresponding
#' to 1 ms).
#' @param granularity Numeric. Mean energy
#' autocorrelation is calculated for intervals of this duration (in ms).
#' @param param Numeric. Voicing onset
#' is predicted by comparing the mean autocorrelation in intervals of size
#' `granularity` after the predicted closure
#' to the interval with the highest autocorrelation,
#' and finding the first interval where mean autocorrelation is at this
#' proportion of the most autocorrelated interval.
#'
#' @return A numeric value giving the sample number of predicted voicing onset.
#' @seealso This function is typically called by either [positiveVOT] or
#' [negativeVOT] with `vo_method='acf'`.
#' These are in turn typically called when using
#' [VOT2newTG] or [addVOT2TG] to add predicted VOT to a Praat TextGrid or to an
#' EMU database with [addVOT2emuDB].
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vd', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#' step <- 0.001*sr
#'
#' clo <- find_closure(sig, sr)
#' vo <- vo_acf(sig, clo, length(sig), clo, step, 1.2, 0.9)
vo_acf <- function(sound, start, end, baseline, step, granularity, param) {
  sq_vo <- seq(from=start, to=end, by=(step*granularity))
  mu_acf <- c()
  i <- 1
  for (s in sq_vo) {
    acf <- stats::acf(sound[s:(s+(step*granularity)-1)], plot=F,
                      na.action=stats::na.pass)
    mu_acf[i] <- mean(acf$acf)
    i <- i+1
  }

  hi_acf <- which(mu_acf > max(mu_acf, na.rm=T)*param)
  if (is.na(hi_acf[2])) {
    vo <- (baseline + (hi_acf[1]*(step*granularity)) + step)
  } else {
    vo <- (baseline + (hi_acf[2]*(step*granularity)) + step)
  }

  return(vo)
}
