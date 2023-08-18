#' Predict negative voice onset time and plot wave with segmented pre-voicing
#'
#' Negative voice onset time and the location of pre-voicing is predicted for a
#' sound snippet containing a pre-voiced stop. Optionally, the wave is plotted
#' showing the predicted location of pre-voicing.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param sr Integer; sample rate of `sound`.
#' @param vo_method String giving the method for predicting the onset of voicing.
#' There are two legal options:
#' * `acf` (default). The onset of voicing is predicted from the mean degree
#' of energy autocorrelation in short intervals of `sound` relative to the most
#' autocorrelated interval during the vowel (see `vo_param`). If audio quality
#' is sufficiently high and stable, `acf` usually gives the best results.
#' * `f0`. The onset of voicing is predicted using the pitch tracking
#' algorithm implemented in [phonTools::pitchtrack].
#' @param closure_interval Numeric. The rough location of the voicing onset is
#' predicted by first finding the quietest interval in the first half of `sound`;
#' `closure_interval` determines the duration of intervals to look for in ms.
#' Default is `10`.
#' @param vo_granularity Numeric, only used when `vo_method='acf'`. Mean energy
#' autocorrelation is calculated for intervals of this duration (in ms). Default
#' is `1.2`.
#' @param vo_param Numeric, only used when `vo_method='acf'`. Voicing onset
#' is predicted by comparing the mean autocorrelation in intervals of size
#' `vo_granularity` after the predicted closure
#' to the interval with the highest autocorrelation,
#' and finding the first interval where mean autocorrelation is at this
#' proportion of the most autocorrelated interval. Default is `0.9`.
#' @param f0_wl Numeric, only used when `vo_method='f0'`. The length of the
#' analysis window in ms passed on to [phonTools::pitchtrack]. Default is `50`.
#' @param f0_minacf Numeric, only used when `vo_method='f0'`. Autocorrelation
#' values below this are ignored by [phonTools::pitchtrack]. Default is `0.5`.
#' @param vo_only Logical; if `TRUE`, only voicing onset is predicted.
#' @param rel_method String giving the method for predicting the location of the
#' burst. There are two legal options:
#' * `transient` (default). The location of the burst is predicted by searching
#' for the transient phase. [phonTools::spectralslice] is used to generate
#' FFT spectra of short portions of `sound` after the predicted onset of voicing,
#' and the smoothest spectrum is predicted to be the location of the transient.
#' If audio quality is sufficiently high and stable, `transient` usually gives
#' the best results.
#' * `amplitude`. The location of the burst is predicted by searching for a
#' sudden increase in amplitude (which usually signifies the beginning of the
#' vowel, which in pre-voiced stops should be very close to the burst).
#' In practice, this is done by recording the highest amplitude in up to 500
#' 1ms chunks of `sound` after the predicted voicing onset, smoothing the
#' resulting time series using a discrete cosine transformation using
#' [emuR::dct], and locating the portion of the resulting function with the
#' highest velocity.
#' @param plot Logical; should the results be plotted? Default is `TRUE`.
#' @param params_list A named list containing the above parameters for
#' determining the onset of voicing. Usually returned by [neg_setParams()], but
#' can also be written by hand. Should have the following named elements (order
#' not important; see details above):
#' * `closure_interval`
#' * `vo_granularity`
#' * `vo_param`
#' * `vo_method`
#' * `f0_wl`
#' * `f0_minacf`
#' * `rel_method`
#' If a `params_list` is provided, it will override any corresponding
#' parameters in the function call.
#'
#' @return A named list with the following elements:
#' * `vo` Gives the sample number where the onset of voicing is predicted.
#' * `rel` Gives the sample number where the stop burst is predicted.
#' * `vot` Gives the predicted voice onset time in ms.
#' * `voi_int` Gives a number with the predicted duration of voicing in
#' `sound` in ms. Potentially used by the [getVOT] function to predict whether
#' voice onset time is positive or negative.
#' @seealso This function is called by the more general and somewhat less
#' flexible function [getVOT], which is in turn called by the functions
#' [VOT2newTG], [addVOT2TG], and [addVOT2emuDB] for bulk analyzing sounds
#' and annotating the output in new or existing TextGrids (for users of the
#' acoustic analysis software Praat) or an EMU database (for users of [emuR]).
#' The parameters can be be somewhat obscure and the default parameters do not
#' necessarily scale particularly well. If a small annotated data set
#' is provided, [neg_setParams] can be used to tune the parameters to your data
#' and annotation style.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vd', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#'
#' #Default settings do not give very good results for this token
#' negativeVOT(sound=sig, sr=sr)
#'
#' #Here with some optimized parameters
#' negativeVOT(sound=sig, sr=sr, closure_interval=3, vo_granularity=1.55,
#' vo_param=0.73)
negativeVOT <- function(sound, sr,
                        vo_method='acf',
                        closure_interval = 10,
                        vo_granularity = 1.2,
                        vo_param=0.9,
                        f0_wl=50, f0_minacf=0.5,
                        vo_only=FALSE,
                        rel_method='transient',
                        plot=TRUE,
                        params_list=NULL) {

  if (length(params_list) > 0) {
    named_params <- names(params_list)
    if ('closure_interval' %in% named_params) {
      closure_interval <- params_list[['closure_interval']]
    }
    if ('rel_method' %in% named_params) {
      release_param <- params_list[['release_param']]
    }
    if ('vo_granularity' %in% named_params) {
      vo_granularity <- params_list[['vo_granularity']]
    }
    if ('vo_param' %in% named_params) {
      vo_param <- params_list[['vo_param']]
    }
    if ('vo_method' %in% named_params) {
      vo_method <- params_list[['vo_method']]
    }
    if ('f0_wl' %in% named_params) {
      f0_wl <- params_list[['f0_wl']]
    }
    if ('f0_minacf' %in% named_params) {
      f0_minacf <- params_list[['f0_minacf']]
    }
    if ('rel_method' %in% named_params) {
      rel_method <- params_list[['rel_method']]
    }
  }

  step <- round(1*(sr/1000))

  if (!(vo_method %in% c('acf', 'f0'))) {
    stop('vo_method must be either acf or f0')
  }
  if (!(rel_method %in% c('transient', 'amplitude'))) {
    stop('rel_method must be either transient or amplitude')
  }

  if (vo_method == 'acf') {
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
    clo <- (which(mean_amp==min(mean_amp))[1] - 0.5) * sr * ci

    sq_vo <- seq(from=clo, to=length(sound), by=(step*vo_granularity))
    mu_acf <- c()
    i <- 1
    for (s in sq_vo) {
      acf <- stats::acf(sound[s:(s+(step*vo_granularity)-1)], plot=F,
                        na.action=stats::na.pass)
      mu_acf[i] <- mean(acf$acf)
      i <- i+1
    }

    hi_acf <- which(mu_acf > max(mu_acf, na.rm=T)*vo_param)
    if (is.na(hi_acf[2])) {
      f0_start <- (clo + (hi_acf[1]*(step*vo_granularity)) + step)
    } else {
      f0_start <- (clo + (hi_acf[2]*(step*vo_granularity)) + step)
    }

  }

  f0 <- phonTools::pitchtrack(sound, fs=sr, show=F, minacf=f0_minacf,
                              windowlength=f0_wl)

  if (vo_method=='f0') {

    if (length(f0$time) > 0) {
      if ((f0$time[2] - f0$time[1]) > 5) {
        f0_start <- f0$time[2] * (sr/1000)
      } else {
        f0_start <- f0$time[1] * (sr/1000)
      }
    } else {f0_start <- NA}

  }

  if (vo_only) {
    if (plot) {
      plot(sound, type='l', x=seq(0, length(sound)/sr, length.out=length(sound)),
           xlab='Time (s)',
           ylab='Amplitude')
      graphics::abline(v=f0_start/sr, col='red')
    }
    return(list(
      vo = f0_start,
      rel = NA,
      vot = NA,
      voi_int = NA
    ))
  } else {
    if (rel_method=='transient') {
      srch <- sound[f0_start:(length(sound)/2)]
      smoothness <- c()

      for (i in 10:(round(length(srch)/step)*10)) {
        stp <- f0_start + ((i*step) / 10)
        slice_this <- sound[(stp-step+1):stp]
        if (sr > 16000) {
          slice_this <- seewave::resamp(slice_this, f=sr, g=16000)[,1]
        }
        slice <- phonTools::spectralslice(slice_this, show=F)[,2]
        smoothness[i-9] <- stats::sd(diff(slice))
      }

      min_smooth <- which(smoothness == min(smoothness))
      pred_rel <- f0_start + (((min_smooth+9)*step) / 10)

    } else if (rel_method=='amplitude') {
      sq_rel <- seq(from=f0_start, to=f0_start+(step*500), by=step)
      max_amp <- c()
      i <- 1
      for (s in sq_rel) {
        max_amp[i] <- max(sound[s:(s+(step)-1)])
        i <- i+1
      }
      if (length(which(is.na(max_amp))) > 0) {
        max_amp <- max_amp[-which(is.na(max_amp))]
      }
      dct_fit <- emuR::dct(max_amp, m=25, fit=T)
      diff_amp <- diff(dct_fit)
      pred_rel <- sq_rel[which.max(diff_amp)] - (step*5)
    }

    if (plot) {
      plot(sound, type='l', x=seq(0, length(sound)/sr, length.out=length(sound)),
           xlab='Time (s)',
           ylab='Amplitude')
      graphics::abline(v=pred_rel/sr, col='red')
      graphics::abline(v=f0_start/sr, col='red')
    }

    return(list(
      vo = f0_start,
      rel = pred_rel,
      vot = -round((pred_rel-f0_start)/sr, 4) * 1000,
      voi_int = length(f0$time) * 2
    ))
  }
}
