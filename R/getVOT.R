#' Predict positive voice onset time and plot wave with segmented stop release
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param sr Integer; sample rate of `sound`.
#' @param plot Logical; if `TRUE` (default), will plot wave with segmented
#' stop release.
#' @param closure_interval Numeric. The rough location of the stop closure is
#' predicted by finding the quietest interval in the first half of `sound`;
#' `closure_interval` determines the duration of intervals to look for.
#' Default is `10`.
#' @param release_param Numeric. The location of the burst is predicted by
#' searching for the first instance after the
#' estimated closure where a 1 ms interval has
#' an energy equal to or higher than this proportion of the maximum energy in
#' `sound`. Default is `15`.
#' @param vo_method String giving the method for predicting the onset of voicing.
#' There are two legal options:
#' * `acf` (default). The beginning of voicing is predicted from the mean degree
#' of energy autocorrelation in short intervals of `sound` relative to the most
#' autocorrelated interval during the vowel. This method seems to work best
#' when `sound` has a high sample rate of e.g. 44.1 kHz or 48 kHz.
#' * `f0`. The beginning of voicing is predicted using the pitch tracking
#' algorithm implemented in `phonTools::pitchtrack()`. This method seeks to work
#' best when `sound` has a low sample rate of e.g. 16 kHz.
#' @param vo_granularity Numeric, only used when `vo_method='acf'`. Mean energy
#' autocorrelation is calculated for intervals of this duration (in ms). Default
#' is `1`.
#' @param vo_param Numeric, only used when `vo_method='acf'`. Voicing onset
#' is predicted by comparing the mean autocorrelation in intervals throughout
#' the release and vowel to the interval with the highest autocorrelation,
#' and finding the first interval where mean autocorrelation is at this
#' proportion of the most autocorrelated interval. Default is `0.85`.
#' @param f0_wl Numeric, only used when `vo_method='f0'`. The length of the
#' analysis window in ms passed on to `phonTools::pitchtrack()`. Default is `30`,
#' much flower than the `phonTools::pitchtrack()` default of `50`, which seems
#' to give better less conservative results when searching for the first pulses
#' after a stop burst.
#' @param f0_minacf Numeric, only used when `vo_method='f0'`. Autocorrelation
#' values below this are ignored by `phonTools::pitchtrack()`. Default is `0.5`.
#'
#' @return
#' @export
#'
#' @examples
getVOT <- function(sound, sr, plot=TRUE,
                   closure_interval = 10,
                   release_param = 15,
                   vo_method = 'acf',
                   vo_granularity = 1,
                   vo_param = 0.85,
                   f0_wl = 30, f0_minacf = 0.5) {

  ci <- closure_interval/1000

  dur <- length(sound)/sr
  sqlen <- ceiling(length(sound)/2 / sr / ci)
  sq_clo <- seq(from=sr*ci, to=sqlen*sr*ci, by=sr*ci)

  mean_amp <- c()
  i <- 1
  for (s in sq_clo) {
    mean_amp[i] <- mean(abs(sound[(s-(sr*ci)+1):s]))
    i <- i+1
  }
  clo <- (which(mean_amp==min(mean_amp)) - 0.5) * sr * ci

  step <- 0.001*sr
  sq_rel <- seq(from=clo[1], to=length(sound), by=step)
  max_amp <- c()
  i <- 1
  for (s in sq_rel) {
    max_amp[i] <- max(sound[s:(s+step-1)])
    i <- i+1
  }

  spike_size <- max(sound)/release_param
  spike <- which(max_amp > spike_size)
  rel <- (clo[1] + ((spike[1])*step))-step

  if (is.na(rel)) {
    return(list(
      rel = clo,
      vo = clo+(0.02*sr),
      vot = 'NA'
    ))
  } else {
    if (!(vo_method %in% c('acf', 'f0'))) {
      stop('Legal parameters for vo_method are acf and f0.')
    }

    if (vo_method == 'acf') {
      sq_vo <- seq(from=rel+(step*5), to=rel+(step*200), by=(step*vo_granularity))
      mu_acf <- c()
      i <- 1
      for (s in sq_vo) {
        acf <- stats::acf(sound[s:(s+(step*vo_granularity)-1)], plot=F, na.action=na.pass)
        mu_acf[i] <- mean(acf$acf)
        i <- i+1
      }

      hi_acf <- which(mu_acf > max(mu_acf, na.rm=T)*vo_param)
      if (is.na(hi_acf[2])) {
        vo <- (rel + (hi_acf[1]*(step*vo_granularity)) + step)
      } else {
        vo <- (rel + (hi_acf[2]*(step*vo_granularity)) + step)
      }


    } else if (vo_method == 'f0') {
      f0 <- phonTools::pitchtrack(sound[rel:length(sound)], fs=sr, show=F,
                                  windowlength=f0_wl, minacf=f0_minacf)
      vo <- round(f0$time[1]) * step + rel

      if (is.na(vo)) {
        return(list(
          rel = rel,
          vo = rel+(0.02*sr),
          vot = 'NA'
        ))
      }
    }

    vot <- round((vo-rel)/sr, 4) * 1000

    if(plot){
      plot(y=sound, x=seq(0, dur, length.out=length(sound)), type='l',
           xlab='Time (s)',
           ylab='Amplitude')
      abline(v=rel/sr, col='red')
      abline(v=vo/sr, col='red')
    }

    return(list(
      rel = rel,
      vo = vo,
      vot = vot
    ))
  }

}
