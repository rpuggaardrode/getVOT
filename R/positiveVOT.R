#' Predict positive voice onset time and plot wave with segmented stop release
#'
#' Positive voice onset time and the location of pre-voicing is predicted for a
#' sound snippet containing a voiceless stop. Optionally, the wave is plotted
#' showing the predicted location of the stop release.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param sr Integer; sample rate of `sound`.
#' @param closure_interval Numeric. The rough location of the stop burst is
#' predicted by first finding the quietest interval in the first half of `sound`;
#' `closure_interval` determines the duration of intervals to look for in ms.
#' Default is `10`.
#' @param release_param Numeric. The location of the burst is predicted by
#' searching for the first instance after the
#' estimated closure where a 1 ms interval has
#' an energy equal to or higher than this proportion of the maximum energy in
#' `sound`. Default is `15`.
#' @param vo_method String giving the method for predicting the onset of voicing.
#' There are two legal options:
#' * `acf` (default). The onset of voicing is predicted from the mean degree
#' of energy autocorrelation in short intervals of `sound` relative to the most
#' autocorrelated interval during the vowel (see `vo_param`). If audio quality
#' is sufficiently high and stable, `acf` usually gives the best results.
#' * `f0`. The onset of voicing is predicted using the pitch tracking
#' algorithm implemented in [phonTools::pitchtrack].
#' * `zcr`. The onset of voicing is predicted from the zero crossing rate.
#' ZCR is calculated in
#' 5 ms windows with 80% overlap, and the resulting time series is smoothed
#' using a discrete cosine transformation with the number of coefficients set to
#' 1/10 of the number of windows. The first ZCR value below the threshold set by
#' `zcr_min` is predicted to correspond to the voicing onset, unless there is
#' an adjacent gap in low ZCR values.
#' @param vo_granularity Numeric, only used when `vo_method='acf'`. Mean energy
#' autocorrelation is calculated for intervals of this duration (in ms). Default
#' is `1`.
#' @param vo_param Numeric, only used when `vo_method='acf'`. Voicing onset
#' is predicted by comparing the mean autocorrelation in intervals of size
#' `vo_granularity` after the predicted closure
#' to the interval with the highest autocorrelation,
#' and finding the first interval where mean autocorrelation is at this
#' proportion of the most autocorrelated interval. Default is `0.85`.
#' @param f0_wl Numeric, only used when `vo_method='f0'`. The length of the
#' analysis window in ms passed on to [phonTools::pitchtrack]. Default is `30`,
#' much lower than the [phonTools::pitchtrack] default of `50`, which seems
#' to give better, less conservative results when searching for the first pulses
#' after a stop burst.
#' @param f0_minacf Numeric, only used when `vo_method='f0'`. Autocorrelation
#' values below this are ignored by [phonTools::pitchtrack]. Default is `0.5`.
#' @param burst_only Logical; if `TRUE`, only burst location is predicted.
#' @param vo_only Boolean; default is `FALSE`. Can be set to `TRUE` if
#' `sign='positive'`, and the data is already aligned such that
#' intervals are aligned to the release. In this case, the burst location is
#' set as the beginning of the interval, and only voicing onset location is
#' predicted.
#' @param rel_offset Numeric, default is `0`. If `vo_only=TRUE`, the algorithm
#' may perform poorly if there is periodicity from e.g. voicing bleed early on
#' in the interval.
#' `rel_offset` tells `getVOT` how much of the initial portion of an interval
#' to ignore when looking for voicing onset (in seconds).
#' @param f0_first Logical; if `TRUE`, the function does not start out by
#' estimating the location of the closure, but rather uses [phonTools::pitchtrack]
#' to find the longest period of consecutive voicing, and searches for the burst
#' in the 200ms chunk prior to the onset of that period.
#' @param zcr_min Numeric, only used when `vo_method='zcr'`. Voicing is
#' predicted to begin when ZCR values after the predicted release
#' are consistently below this threshold.
#' Default is `0.02`.
#' @param plot Logical; should the results be plotted? Default is `TRUE`.
#' @param params_list A named list containing the above parameters for
#' determining the onset of voicing. Usually returned by [pos_setParams()], but
#' can also be written by hand. Should have the following named elements (order
#' not important; see details above):
#' * `closure_interval`
#' * `vo_granularity`
#' * `vo_param`
#' * `vo_method`
#' * `f0_wl`
#' * `f0_minacf`
#' * `release_param`
#' * `f0_first`
#' If a `params_list` is provided, it will override any corresponding
#' parameters in the function call.
#'
#' @return A named list with the following elements:
#' * `rel` Gives the sample number where the stop burst is predicted.
#' * `vo` Gives the sample number where the onset of voicing is predicted.
#' * `vot` Gives the predicted voice onset time in ms.
#' * `spike_size` Gives the amplitude of the burst relative to the highest
#' measured amplitude in the signal. Potentially used by the [getVOT] function
#' to predict whether voice onset time is positive or negative.
#' @seealso This function is called by the more general and somewhat less
#' flexible function [getVOT], which is in turn called by the functions
#' [VOT2newTG], [addVOT2TG], and [addVOT2emuDB] for bulk analyzing sounds
#' and annotating the output in new or existing TextGrids (for users of the
#' acoustic analysis software Praat) or an EMU database (for users of [emuR]).
#' The parameters can be be somewhat obscure and the default parameters do not
#' necessarily scale particularly well. If a small annotated data set
#' is provided, [pos_setParams] can be used to tune the parameters to your data
#' and annotation style.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#'
#' positiveVOT(sound=sig, sr=sr)
positiveVOT <- function(sound, sr,
                        closure_interval = 10,
                        release_param = 15,
                        vo_method = 'acf',
                        vo_granularity = 1,
                        vo_param = 0.85,
                        f0_wl = 30, f0_minacf = 0.5,
                        zcr_min = 0.05,
                        burst_only=FALSE,
                        vo_only=FALSE,
                        rel_offset=0,
                        f0_first=FALSE,
                        plot=TRUE,
                        params_list=NULL) {

  if (length(params_list) > 0) {
    named_params <- names(params_list)
    if ('closure_interval' %in% named_params) {
      closure_interval <- params_list[['closure_interval']]
    }
    if ('release_param' %in% named_params) {
      release_param <- params_list[['release_param']]
    }
    if ('vo_granularity' %in% named_params) {
      vo_granularity <- params_list[['vo_granularity']]
    }
    if ('vo_param' %in% named_params) {
      vo_param <- params_list[['vo_param']]
    }
    if ('f0_first' %in% named_params) {
      f0_first <- params_list[['f0_first']]
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
    if ('zcr_min' %in% named_params) {
      zcr_min <- params_list[['zcr_min']]
    }
  }

  if (!(vo_method %in% c('acf', 'f0', 'zcr'))) {
    stop('vo_method must be either acf or f0')
  }

  ci <- closure_interval/1000
  step <- 0.001*sr
  dur <- length(sound)/sr

  if (f0_first & !vo_only) {
    f0 <- phonTools::pitchtrack(sound, fs=sr, show=F)
    time_vec <- f0$time

    adj_f0 <- c()
    f0_seq <- seq(0, dur*1000, 100)
    for (s in 1:length(f0_seq)) {
      adj_f0[s] <- length(which(time_vec > f0_seq[s]+1 & time_vec < f0_seq[s]+100))
    }
    if (max(adj_f0) > 40) {
      many_adj_f0 <- (which(adj_f0 > 40)[1] * 100) - 50
    } else {
      many_adj_f0 <- (which.max(adj_f0) * 100) - 50
    }

    v_loc <- which.min(abs(many_adj_f0 - time_vec))
    diff_vec <- as.numeric(diff(time_vec) <= 2)
    vo_index <- rev(which(diff_vec[1:v_loc] == 0))[1] + 1
    if (is.na(vo_index)) {
      vo <- time_vec[1] * step
    } else {
      vo <- time_vec[vo_index] * step
    }

    if (vo-(sr*0.2) > 0) {
      sq_rel <- seq(from=vo-(sr*0.2), to=vo, by=step)
    } else {
      sq_rel <- seq(from=vo-(sr*0.05), to=vo, by=step)
    }
    max_amp <- c()
    i <- 1
    for (s in sq_rel) {
      max_amp[i] <- max(abs(sound[s:(s+step-1)]))
      i <- i+1
    }

    spike_size <- max(abs(sound))/release_param
    spike <- which(max_amp > spike_size)
    rel <- (vo-(step*200) + ((spike[1])*step))-step
  } else if (!f0_first & !vo_only) {
    sqlen <- ceiling(length(sound)/2 / sr / ci)
    sq_clo <- seq(from=sr*ci, to=sqlen*sr*ci, by=sr*ci)

    mean_amp <- c()
    i <- 1
    for (s in sq_clo) {
      mean_amp[i] <- mean(abs(sound[(s-(sr*ci)+1):s]))
      i <- i+1
    }
    clo <- (which(mean_amp==min(mean_amp)) - 0.5) * sr * ci

    sq_rel <- seq(from=clo[1], to=length(sound), by=step)
    max_amp <- c()
    i <- 1
    for (s in sq_rel) {
      max_amp[i] <- max(abs(sound[s:(s+step-1)]))
      i <- i+1
    }

    spike_size <- max(abs(sound))/release_param
    spike <- which(max_amp > spike_size)
    rel <- (clo[1] + ((spike[1])*step))-step
  } else if (vo_only) {
    rel <- 0 + (rel_offset*sr)
    spike_size <- NA
  }

  if (!vo_only) {
    if (is.na(rel)) {
      if (f0_first) {
        return(list(
          rel = vo-(sr*0.02),
          vo = vo,
          vot = 'NA'
        ))
      } else {
        return(list(
          rel = clo,
          vo = clo+(0.02*sr),
          vot = 'NA'
        ))
      }
    }
  }


  if (burst_only) {
    if (plot) {
      plot(sound, type='l', x=seq(0, length(sound)/sr, length.out=length(sound)),
           xlab='Time (s)',
           ylab='Amplitude')
      graphics::abline(v=rel/sr, col='red')
    }
    return(list(
      rel = rel,
      vo = NA,
      vot = NA,
      spike=spike_size
    ))
  } else {
    if (vo_method == 'acf') {
      sq_vo <- seq(from=rel+(step*5), to=rel+(step*200), by=(step*vo_granularity))
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

    } else if (vo_method == 'zcr') {
      vo_srch <- sound[rel:length(sound)]
      zcr_ts <- seewave::zcr(vo_srch, sr, wl=sr*0.005, ovlp=80, plot=F)
      zcr_smooth <- emuR::dct(zcr_ts[,'zcr'], m=length(zcr_ts)/20, fit=T)
      zcr_smooth[1:10] <- NA

      if (min(zcr_smooth, na.rm=T) < zcr_min) {
        low_zcr <- which(zcr_smooth < zcr_min)
        if (max(diff(low_zcr)) < 40) {
          f0_start <- (zcr_ts[low_zcr[1], 'time'] * sr)
        } else {
          first_voi <- which(diff(low_zcr) > 40)[1] + 1
          if (first_voi > 20 | is.na(first_voi)) {
            f0_start <- (zcr_ts[low_zcr[1], 'time'] * sr)
          } else {
            f0_start <- (zcr_ts[low_zcr[first_voi], 'time'] * sr)
          }
        }
      } else {
        f0_start <- (min(zcr_ts[,'zcr']) * sr)
      }

      vo <- f0_start + rel
      names(vo) <- NULL

    }

    if (rel_offset > 0) rel <- rel - (rel_offset*sr)
    vot <- round((vo-rel)/sr, 4) * 1000

    if(plot){
      plot(y=sound, x=seq(0, dur, length.out=length(sound)), type='l',
           xlab='Time (s)',
           ylab='Amplitude')
      graphics::abline(v=rel/sr, col='red')
      graphics::abline(v=vo/sr, col='red')
    }

    return(list(
      rel = rel,
      vo = vo,
      vot = vot,
      spike = spike_size
    ))
  }

}
