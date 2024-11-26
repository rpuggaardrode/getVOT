#' Predict voicing onset in a stop based on strength of harmonic excitation
#'
#' Predict the location of voicing in a short sound file by passing a
#' differenced, downsampled version of the waveform through a cascade of zero
#' frequency filters to derive the strength of harmonic excitation (SoE), and
#' choosing the first zero crossing where the SoE is above a certain threshold.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param start Numeric value giving the first sample of the search space in
#' the sound signal.
#' @param sr Integer; sample rate of `sound`.
#' @param min_soe Numeric, only used when `vo_method='soe'`. Voicing is
#' predicted to begin by the first zero crossing with this SoE value in a
#' filtered version of the waveform. Default is `0.01`.
#'
#' @return A numeric value giving the sample number of predicted voicing onset.
#' @seealso This function is typically called by either [positiveVOT] or
#' [negativeVOT] with `vo_method='soe'`.
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
#' vo <- vo_soe(sig, clo, sr)
vo_soe <- function(sound, start, sr, min_soe = 0.01) {

  vo_srch <- sound[start:length(sound)]

  waveObj <- tuneR::Wave(sound, samp.rate=sr, bit=16)
  new_sr <- 16000
  sound <- tuneR::downsample(waveObj, new_sr)@left
  wid <- 201

  sig <- diff(sound, 1)
  len <- length(sig)

  zfr1_filt <- signal::filter(1, c(1, -2*0.999, 0.999^2), sig)
  a <- signal::filter(rep(1, wid)/wid, 1, zfr1_filt)
  abegin <- cumsum(sig[1:(wid-2)])
  abegin <- abegin[seq(1, length(abegin), by=2)] / seq(1, wid-2, by=2)
  aend <- cumsum(sig[len:(len-wid+3)])
  aend <- aend[seq(length(aend), 1, by=-2)] / seq(wid-2, 1, by=-2)
  a <- c(abegin, a[wid:length(a)], aend)
  zfr1_trendRem <- zfr1_filt - a

  zfr2_filt <- signal::filter(1, c(1, -2*0.999, 0.999^2), zfr1_trendRem)
  a <- signal::filter(rep(1, wid)/wid, 1, zfr2_filt)
  abegin <- cumsum(zfr2_filt[1:(wid-2)])
  abegin <- abegin[seq(1, length(abegin), by=2)] / seq(1, wid-2, by=2)
  aend <- cumsum(zfr2_filt[len:(len-wid+3)])
  aend <- aend[seq(length(aend), 1, by=-2)] / seq(wid-2, 1, by=-2)
  a <- c(abegin, a[wid:length(a)], aend)

  zfr_out <- zfr2_filt - a

  z <- 0.95*zfr_out[1:(length(zfr_out)-wid)]/
    max(abs(zfr_out[1:(length(zfr_out)-wid)]))
  z1 <- c(NA, z[1:length(z)-1])
  tf <- z1 > 0 & z<=0
  pulses <- which(tf)
  dx <- c()
  for (i in 1:length(pulses)) {
    dx[i] <- stats::lm(z1[(pulses[i]-1):(pulses[i]+1)]~c(1:3))$coefficients[2]
  }
  soe <- -dx[-1]
  if (any(soe > min_soe)) {
    vo_new_sr <- pulses[which(soe > min_soe)][1]
  } else {
    vo_new_sr <- pulses[which.max(soe)]
  }

  vo <- (vo_new_sr / new_sr) * sr
  vo <- start + vo
  return(vo)

}
