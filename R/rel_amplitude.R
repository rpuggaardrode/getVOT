#' Predict location of stop release based on relative amplitude
#'
#' Predicts the location of a stop release in a short chunked sound snippet
#' or differenced speech signal
#' by searching for which chunk is the first to have an amplitude peak above
#' a certain threshold relative to the highest amplitude in the snippet.
#' Optionally, the peaks are treated as a smoothed time series, in which case
#' the release is assumed to occur at the velocity peak of the smooth function.
#' The latter option is mainly used for estimating releases in pre-voiced stops.
#'
#' @param sound Numeric vector corresponding to a sound wave or differenced
#' speech signal.
#' @param start Numeric value giving the first sample of the search space in
#' the sound signal.
#' @param sq_rel Numeric vector determining the locations of chunks in the
#' sound snippet.
#' @param step The duration of each chunk in terms of samples (corresponding
#' to 1 ms).
#' @param release_param Numeric. The location of the burst is predicted by
#' searching for the first instance after the
#' estimated closure where a 1 ms interval has
#' an energy equal to or higher than this proportion of the maximum energy in
#' `sound`. Default is `15`.
#' @param smooth Boolean; should the release be estimated by searching for the
#' velocity peak in a smoothed time series of amplitude peaks? Default is
#' `FALSE`.
#'
#' @return When `smooth=FALSE`, a named list containing `rel`, the sample number
#' where the release is predicted to be, and `spike_size`, the size of the
#' amplitude peak at the predicted release. When `smooth=TRUE`, returns only
#' the sample number where the release is predicted to be.
#' @seealso This function is typically called by either [positiveVOT] or
#' [negativeVOT] with `rel_method='amplitude'`, `rel_method='raw'`, or
#' `rel_method='diff'`. These are in turn typically called when using
#' [VOT2newTG] or [addVOT2TG] to add predicted VOT to a Praat TextGrid or to an
#' EMU database with [addVOT2emuDB].
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#' step <- 0.001*sr
#'
#' clo <- find_closure(sig, sr)
#' sq_rel <- seq(from=clo, to=length(sig), by=step)
#' rel <- rel_amplitude(sig, clo, sq_rel, step)
rel_amplitude <- function(sound, start, sq_rel, step,
                          release_param = 15,
                          smooth = FALSE) {
  max_amp <- c()
  i <- 1
  for (s in sq_rel) {
    max_amp[i] <- max(abs(sound[s:(s+step-1)]))
    i <- i+1
  }
  if (length(which(is.na(max_amp))) > 0) {
    max_amp <- max_amp[-which(is.na(max_amp))]
  }

  if (smooth) {
    dct_fit <- emuR::dct(max_amp, m=25, fit=T)
    diff_amp <- diff(dct_fit)
    pred_rel <- sq_rel[which.max(diff_amp)] - (step*5)
    return(pred_rel)
  } else {
    spike_size <- max(abs(sound))/release_param
    spike <- which(max_amp > spike_size)
    pred_rel <- (start + ((spike[1])*step))-step
    return(list(rel = pred_rel, spike_size = spike_size))
  }
}
