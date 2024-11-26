#' Optimize parameters for predicting negative voice onset time
#'
#' If the user provides a few examples of WAV/TextGrid pairs with hand-annotated
#' voice onset time of stops, this function will attempt to predict VOT with
#' a large range of parameters, and select those parameters that on average
#' minimize prediction error the most. These parameters can then be used for
#' the rest of the data set.
#'
#' @param directory A string giving the location of a directory with WAV/TextGrid
#' pairs. Ideally, each WAV file should contain just one syllable with a stop,
#' where the stop closure is placed somewhere in the first half of the sound file.
#' The TextGrid tier containing hand-segmented VOT should be named `vot`.
#' Default is the current working directory.
#' @param when_avoid_soe Numeric; default is `10`. Per default, voicing onset is
#' predicted with [negativeVOT] using `method='soe'` (see [vo_soe]). This is the
#' preferred method because it is much faster than the alternative methods.
#' For this reason, other methods are only tested if the average prediction
#' error with `method='acf'` is above a certain threshold. `when_avoid_soe` sets
#' this threshold ,i.e. by default, other methods are tested when the prediction
#' error with `method='soe'` is above 10 ms.
#' @param when_use_f0 Numeric; default is `10`. Per default, voicing onset is
#' predicted with [negativeVOT] using `method='soe'`, alternatively
#' `method='acf'` or `method='zcr'`. These methods are preferred to `method='f0'`
#' because it is quite slow.  For
#' this reason, parameters with `method='f0'` are only tested if the average
#' prediction error with other methods are above a certain threshold.
#' `when_use_f0` sets this threshold; i.e. by default, parameters for
#' `method='f0'` are tested when the average prediction error with the best
#' alternative method is above 10 ms.
#' @param verbose Logical; should messages be printed in the console indicating
#' the progress? Default is `TRUE`.
#' @param plot Logical; should predicted VOT using the resulting parameters
#' be plotted against the hand-annotated VOT?
#' @param ... Additional graphical parameters to be passed on to
#' [plot_training_diff].
#'
#' @return A named list with optimized parameters for the training data; see
#' [negativeVOT] for more information.
#' @seealso The returned list can be passed to the [negativeVOT] argument
#' `params_list`, or the argument `neg_params_list` in other functions
#' [getVOT], [VOT2newTG], [addVOT2TG], [addVOT2emuDB]. [pos_setParams] is a
#' sister function for optimizing parameters for predicting positive voice
#' onset time. [plot_training_diff] can be used to illustrate how closely the
#' optimized parameters predict the hand-annotated training data.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vd', package='getVOT')
#' neg_setParams(directory=datapath)
neg_setParams <- function(directory='.',
                          when_avoid_soe=10,
                          when_use_f0=10,
                          verbose=TRUE,
                          plot=TRUE,
                          ...) {

  tg_list <- list.files(directory, '*.TextGrid')
  wav_list <- list.files(directory, '*.wav')
  if (!(identical(unlist(strsplit(wav_list, '*.wav')),
                  unlist(strsplit(tg_list, '*.TextGrid'))))) {
    print(paste0('There is not a one-to-one correspondence between sound files ',
                 'and TextGrids in the directory'))
  }

  pred <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Testing voicing onset parameter settings for ', wav_list[snd]))
    }

    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

    i <- 1

    for (ci in c(5, 10, 15, 20)) {
      for (soe in c(0.003, 0.008, 0.013, 0.018)) {
        tmp <- negativeVOT(sig[,1], fs, vo_method='soe', closure_interval=ci,
                           soe_min=soe, vo_only=T, plot=F)

        pred$ci[i] <- ci
        pred$gran[i] <- NA
        pred$p[i] <- NA
        pred$method[i] <- 'soe'
        pred$zcrmin[i] <- NA
        pred$soemin[i] <- soe
        pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
        i <- i+1
      }
    }
  }

  pred <- as.data.frame(pred)
  pred$diff_mean <- rowMeans(pred[7:ncol(pred)])

  if (verbose) {
    print(paste0('On average, the selected parameters for voicing onset after a first pass ',
                 'agree with the training data within a margin of ',
                 round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
  }

  winner <- which.min(pred$diff_mean)
  good_ci <- pred$ci[winner]
  good_gran <- pred$gran[winner]
  good_p <- pred$p[winner]
  good_z <- pred$zcrmin[winner]
  good_soe <- pred$soemin[winner]
  good_method <- pred$method[winner]

  finetune <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Finetuning voicing onset parameter settings for ', wav_list[snd]))
    }

    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

    i <- 1

    for (ci in c(good_ci-2, good_ci-1, good_ci, good_ci+1, good_ci+2)) {
      for (soe in c(good_soe-0.002, good_soe-0.001, good_soe,
                    good_soe+0.001, good_soe+0.002)) {
        tmp <- negativeVOT(sig[,1], fs, vo_method='soe', closure_interval = ci,
                           soe_min = soe, vo_only = T, plot = F)

        finetune$ci[i] <- ci
        finetune$gran[i] <- NA
        finetune$p[i] <- NA
        finetune$method[i] <- 'soe'
        finetune$zcrmin[i] <- NA
        finetune$soemin[i] <- soe
        finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
        i <- i+1
      }
    }
  }

  finetune <- as.data.frame(finetune)
  finetune$diff_mean <- rowMeans(finetune[7:ncol(finetune)])

  if (verbose) {
    print(paste0('On average, the selected voicing onset parameters after finetuning ',
                 'agree with the training data within a margin of ',
                 round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
  }

  winner <- which.min(finetune$diff_mean)
  ci_winner <- finetune$ci[winner]
  gran_winner <- finetune$gran[winner]
  p_winner <- finetune$p[winner]
  z_winner <- finetune$zcrmin[winner]
  soe_winner <- finetune$soemin[winner]
  with_soe <- min(finetune$diff_mean) / (fs/1000)

  if (with_soe > when_avoid_soe) {
    if (verbose) {
      print(paste0('Average agreement with training data below ',
                   when_avoid_soe,
                   ' ms. Trying vo_method=acf and vo_method=zcr'))
    }

    for (snd in 1:length(wav_list)){
      if (verbose) {
        print(paste0('Testing voicing onset parameter settings for ', wav_list[snd]))
      }

      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i <- 1

      for (ci in c(5, 10, 15, 20)) {

        for (gran in c(0.8, 1, 1.2, 1.4, 1.6)) {

          for (p in c(0.75, 0.8, 0.85, 0.9)) {
            tmp <- negativeVOT(sig[,1], fs, vo_method='acf', closure_interval = ci,
                               vo_granularity = gran, vo_param = p, vo_only=T,
                               plot=F)

            pred$ci[i] <- ci
            pred$gran[i] <- gran
            pred$p[i] <- p
            pred$method[i] <- 'acf'
            pred$zcrmin[i] <- NA
            pred$soemin[i] <- NA
            pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
            i <- i+1
          }
        }

        for (z in c(0.02, 0.04, 0.06, 0.08, 0.1)) {
          tmp <- negativeVOT(sig[,1], fs, vo_method='zcr', closure_interval = ci,
                             zcr_min = z, vo_only = T, plot = F)

          pred$ci[i] <- ci
          pred$gran[i] <- NA
          pred$p[i] <- NA
          pred$method[i] <- 'zcr'
          pred$zcrmin[i] <- z
          pred$soemin[i] <- NA
          pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
          i <- i+1
        }
      }

    }

    pred <- as.data.frame(pred)
    pred$diff_mean <- rowMeans(pred[7:ncol(pred)])

    winner <- which.min(pred$diff_mean)
    good_ci <- pred$ci[winner]
    good_gran <- pred$gran[winner]
    good_p <- pred$p[winner]
    good_z <- pred$zcrmin[winner]
    good_soe <- pred$soemin[winner]
    good_method <- pred$method[winner]

    finetune <- list()

    for (snd in 1:length(wav_list)){

      if (verbose) {
        print(paste0('Finetuning voicing onset parameter settings for ', wav_list[snd]))
      }

      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i <- 1

      for (ci in c(good_ci-2, good_ci-1, good_ci, good_ci+1, good_ci+2)) {
        if (good_method == 'acf') {
          for (gran in c(good_gran-0.1, good_gran-0.05, good_gran,
                         good_gran+0.05, good_gran+0.1)) {

            for (p in c(good_p-0.02, good_p-0.01, good_p,
                        good_p+0.01, good_p+0.02)) {
              tmp <- negativeVOT(sig[,1], fs, vo_method='acf', closure_interval = ci,
                                 vo_granularity = gran, vo_param = p, vo_only=T,
                                 plot=F)

              finetune$ci[i] <- ci
              finetune$gran[i] <- gran
              finetune$p[i] <- p
              finetune$method[i] <- 'acf'
              finetune$zcrmin[i] <- NA
              finetune$soemin[i] <- NA
              finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
              i <- i+1
            }
          }
        }

        if (good_method == 'zcr') {
          for (z in c(good_z-0.01, good_z-0.005, good_z, good_z+0.005, good_z+0.01)) {
            tmp <- negativeVOT(sig[,1], fs, vo_method='zcr', closure_interval = ci,
                               zcr_min = z, vo_only = T, plot = F)

            finetune$ci[i] <- ci
            finetune$gran[i] <- NA
            finetune$p[i] <- NA
            finetune$method[i] <- 'zcr'
            finetune$zcrmin[i] <- z
            finetune$soemin[i] <- NA
            finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
            i <- i+1
          }
        }
      }
    }

    finetune <- as.data.frame(finetune)
    finetune$diff_mean <- rowMeans(finetune[6:ncol(finetune)])

    if (verbose) {
      print(paste0('On average, the selected voicing onset parameters after finetuning ',
                   'agree with the training data within a margin of ',
                   round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
    }

    winner <- which.min(finetune$diff_mean)
    ci_winner <- finetune$ci[winner]
    gran_winner <- finetune$gran[winner]
    p_winner <- finetune$p[winner]
    z_winner <- finetune$z[winner]
    soe_winner <- finetune$soe[winner]
    no_f0 <- min(finetune$diff_mean) / (fs/1000)

    if (no_f0 < with_soe) {
      vo_method <- good_method
    } else {
      vo_method <- 'soe'
    }
  } else {
    vo_method <- 'soe'
    no_f0 <- with_soe
  }

  if (no_f0 > when_use_f0) {
    if (verbose) {
      print(paste0('Average agreement with training data below ',
                   when_use_f0,
                   ' ms. Trying vo_method=f0'))
    }
    pred <- list()

    for (snd in 1:length(wav_list)){
      if (verbose) {
        print(paste0('Testing voicing onset parameter settings for ', wav_list[snd]))
      }
      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i = 1

      for (wl in c(30, 40, 50, 60, 70)) {
        for (acf in c(0.4, 0.5, 0.6, 0.7, 0.8)) {
          tmp <- negativeVOT(sig[,1], fs, vo_method='f0',
                             f0_wl = wl, f0_minacf = acf,
                             plot=F, vo_only=T)

          pred$wl[i] <- wl
          pred$acf[i] <- acf
          pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
          i <- i+1
        }
      }
    }

    pred <- as.data.frame(pred)
    pred$diff_mean <- rowMeans(pred[3:ncol(pred)])

    if (verbose) {
      print(paste0('On average, the selected voicing onset parameters after a ',
                   'first pass with method=f0 ',
                   'agree with the training data within a margin of ',
                   round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(pred$diff_mean)
    good_wl <- pred$wl[winner]
    good_acf <- pred$acf[winner]

    finetune <- list()

    for (snd in 1:length(wav_list)){
      if (verbose) {
        print(paste0('Finetuning voicing onset parameter settings with method=f0 for ',
                     wav_list[snd]))
      }
      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i = 1

      for (wl in c(good_wl-5, good_wl-2.5, good_wl, good_wl+2.5, good_wl+5)) {
        for (acf in c(good_acf-0.05, good_acf-0.025, good_acf, good_acf+0.025, good_acf+0.05)) {
          tmp <- negativeVOT(sig[,1], fs, vo_method='f0',
                             f0_wl=wl, f0_minacf=acf,
                             plot=F, vo_only=T)

          finetune$wl[i] <- wl
          finetune$acf[i] <- acf
          finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
          i <- i+1
        }
      }
    }

    finetune <- as.data.frame(finetune)
    finetune$diff_mean <- rowMeans(finetune[3:ncol(finetune)])

    if (verbose) {
      print(paste0('On average, the selected voicing parameters after fin
                   etuning ',
                   'with method=f0 agree with the training data within a margin of ',
                   round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(finetune$diff_mean)[1]
    wl_winner <- finetune$wl[winner]
    acf_winner <- finetune$acf[winner]
    f0_results <- min(finetune$diff_mean) / (fs/1000)

    if (f0_results < no_f0) {
      vo_method <- 'f0'
    } else {
      vo_method <- good_method
    }
  } else {
    vo_method <- vo_method
    wl_winner <- 50
    acf_winner <- 0.5
  }

  pred <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Testing accuracy of different burst detection methods for ',
                   wav_list[snd]))
    }

    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

    tmp <- negativeVOT(sig[,1], fs, vo_method=vo_method, f0_wl=wl_winner,
                              f0_minacf = acf_winner, closure_interval = ci_winner,
                              vo_granularity=gran_winner,
                              vo_param=p_winner, zcr_min=z_winner,
                              rel_method='transient', plot=F)
    pred$rp[1] <- NA
    pred$method[1] <- 'transient'
    pred[[wav_list[snd]]][1] <- abs(tg_t1 - tmp$rel)

    i <- 2

    for (rm in c('amplitude', 'diff')) {
      for (rp in c(9, 12, 15, 18, 21)) {
        tmp <- negativeVOT(sig[,1], fs, vo_method=vo_method, f0_wl=wl_winner,
                           f0_minacf = acf_winner, closure_interval = ci_winner,
                           vo_granularity=gran_winner,
                           vo_param=p_winner, zcr_min=z_winner,
                           rel_method=rm, release_param = rp, plot=F)

        pred$rp[i] <- rp
        pred$method[i] <- rm
        pred[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$rel)
        i <- i+1
      }
    }
  }

  pred <- as.data.frame(pred)
  pred$diff_mean <- rowMeans(pred[3:ncol(pred)])

  if (verbose) {
    print(paste0('On average, the selected burst detection parameters after a ',
                 'first pass ',
                 'agree with the training data within a margin of ',
                 round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
  }
  winner <- which.min(pred$diff_mean)
  good_rp <- pred$rp[winner]
  good_method <- pred$method[winner]

  if (good_method == 'transient') {
    rel_method_winner <- 'transient'
    rp_winner <- NA
  } else {

    finetune <- list()

    for (snd in 1:length(wav_list)){

      if (verbose) {
        print(paste0('Finetuning burst detection parameter settings for ', wav_list[snd]))
      }

      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

      i <- 1

      for (rm in c('amplitude', 'diff')) {
        for (rp in c(good_rp-1, good_rp-0.5, good_rp, good_rp+0.5, good_rp+1)) {
          tmp <- negativeVOT(sig[,1], fs, vo_method=vo_method, f0_wl=wl_winner,
                             f0_minacf = acf_winner, closure_interval = ci_winner,
                             vo_granularity=gran_winner,
                             vo_param=p_winner, zcr_min=z_winner,
                             soe_min = soe_winner,
                             rel_method=rm, release_param = rp, plot=F)

          finetune$rp[i] <- rp
          finetune$method[i] <- rm
          finetune[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$rel)
          i <- i+1
        }
      }
    }

    finetune <- as.data.frame(finetune)
    finetune$diff_mean <- rowMeans(finetune[3:ncol(finetune)])

    if (verbose) {
      print(paste0('On average, the selected burst detection parameters after finetuning ',
                   'agree with the training data within a margin of ',
                   round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(finetune$diff_mean)[1]
    rp_winner <- finetune$rp[winner]
    rel_method_winner <- finetune$method[winner]
  }

  pl <- list(
    closure_interval = ci_winner,
    vo_granularity = gran_winner,
    vo_param = p_winner,
    vo_method = vo_method,
    f0_wl = wl_winner,
    f0_minacf = acf_winner,
    rel_method = rel_method_winner,
    release_param = rp_winner,
    zcr_min = z_winner,
    soe_min = soe_winner
  )

  if (plot) {
    plot_training_diff(directory, pl, ...)
  }

  return(pl)

}
