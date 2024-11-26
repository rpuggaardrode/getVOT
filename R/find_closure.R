#' Predict location of stop closure
#'
#' Predicts the location of a stop closure in a short sound snippet by
#' chunking up the first half of the snippet into equidistant intervals and
#' finding the overall most silent interval.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param sr Integer; sample rate of `sound`.
#' @param closure_interval Numeric. `closure_interval` determines the duration
#' of silent intervals to look for in ms. Default is `10`.
#'
#' @return A numeric value giving the sample number at the centre of the
#' predicted stop closure interval.
#' @seealso This function is typically called by either [positiveVOT] or
#' [negativeVOT], which are in turn typically called when using [VOT2newTG]
#' or [addVOT2TG] to add predicted VOT to a Praat TextGrid or to an
#' EMU database with [addVOT2emuDB].
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#'
#' clo <- find_closure(sig, sr)
find_closure <- function(sound, sr, closure_interval = 10) {

  ci <- closure_interval / 1000

  sqlen <- ceiling(length(sound)/2 / sr / ci)
  sq_clo <- seq(from=sr*ci, to=sqlen*sr*ci, by=sr*ci)

  mean_amp <- c()
  i <- 1
  for (s in sq_clo) {
    mean_amp[i] <- mean(abs(sound[(s-(sr*ci)+1):s]))
    i <- i+1
  }
  clo <- (which(mean_amp==min(mean_amp))[1] - 0.5) * sr * ci

  return(clo)

}
