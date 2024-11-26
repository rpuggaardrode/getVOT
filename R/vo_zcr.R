#' Predict voicing onset in a stop based on zero crossing rate
#'
#' Predict the location of voicing in a short sound file by calculating ZCR in
#' 5 ms windows with 80% overlap, and smoothing the resulting time series
#' using a discrete cosine transformation with the number of coefficients set to
#' 1/10 of the number of windows. The first ZCR value below a certain threshold
#' is predicted to correspond to the voicing onset, unless there is
#' an adjacent gap in low ZCR values.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param start Numeric value giving the first sample of the search space in
#' the sound signal.
#' @param sr Integer; sample rate of `sound`.
#' @param minval Voicing is
#' predicted to begin when ZCR values are consistently below this threshold.
#' @param padNA Boolean; should the first few ZCR values be treated as `NA` in
#' order to avoid predicted voicing onset immediately after the value set by
#' `start`? This is typically necessary when predicting positive VOT using ZCR.
#'
#' @return A numeric value giving the sample number of predicted voicing onset.
#' @seealso This function is typically called by either [positiveVOT] or
#' [negativeVOT] with `vo_method='zcr'`.
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
#'
#' clo <- find_closure(sig, sr)
#' vo <- vo_zcr(sig, clo, sr, 0.02, FALSE)
vo_zcr <- function(sound, start, sr, minval, padNA) {

  vo_srch <- sound[start:length(sound)]
  zcr_ts <- seewave::zcr(vo_srch, sr, wl=sr*0.005, ovlp=80, plot=F)
  zcr_smooth <- emuR::dct(zcr_ts[,'zcr'], m=length(zcr_ts)/20, fit=T)

  if (padNA) zcr_smooth[1:10] <- NA

  if (min(zcr_smooth, na.rm=T) < minval) {
    low_zcr <- which(zcr_smooth < minval)
    if (max(diff(low_zcr)) < 40) {
      vo <- (zcr_ts[low_zcr[1], 'time'] * sr)
    } else {
      first_voi <- which(diff(low_zcr) > 40)[1] + 1
      if (first_voi > 20 | is.na(first_voi)) {
        vo <- (zcr_ts[low_zcr[1], 'time'] * sr)
      } else {
        vo <- (zcr_ts[low_zcr[first_voi], 'time'] * sr)
      }
    }
  } else {
    vo <- (min(zcr_ts[,'zcr']) * sr)
  }

  vo <- vo + start
  names(vo) <- NULL
  return(vo)

}
