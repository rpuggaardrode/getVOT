#' Optimize parameters for predicting positive voice onset time
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
#' @param when_use_f0 Numeric; default is `10`. Per default, voicing onset is
#' predicted with [positiveVOT] using `method='acf'`. This is the preferred method
#' because it is significantly faster than the alternative `method='f0'`. For
#' this reason, parameters with `method='f0'` are only tested if the average
#' prediction error with `method='acf'` is above a certain threshold.
#' `when_use_f0` sets this threshold; i.e. by default, parameters for
#' `method='f0'` are tested when the average prediction error with `method='acf'`
#' is above 10 ms.
#' @param when_use_f0_first Numeric; default is `10`. Per default, the burst is
#' located by first predicting the rough location of the stop closure (i.e.
#' [positiveVOT] is called with `f0_first=FALSE`). This is the preferred method
#' because it is significantly faster than the alternative method, where the
#' burst is predicted by first finding the location of the vowel and moving
#' backwards (`f0_first=TRUE`). For this reason, parameters with `f0_first=TRUE`
#' are only tested if the average prediction error with `f0_first=FALSE` is
#' above a certain threshold. `when_use_f0_first` sets this threshold; by
#' default, parameters with `f0_first=TRUE` are tested when prediction error
#' for half of the training material is above 10 ms.
#' @param verbose Logical; should messages be printed in the console indicating
#' the progress? Default is `TRUE`.
#' @param plot Logical; should predicted VOT using the resulting parameters
#' be plotted against the hand-annotated VOT?
#'
#' @return A named list with optimized parameters for the training data; see
#' [positiveVOT] for more information.
#' @seealso The returned list can be passed to the [positiveVOT] argument
#' `params_list`, or the argument `pos_params_list` in other functions
#' [getVOT], [VOT2newTG], [addVOT2TG], [addVOT2emuDB]. [neg_setParams] is a
#' sister function for optimizing parameters for negative positive voice
#' onset time. [plot_training_diff] can be used to illustrate how closely the
#' optimized parameters predict the hand-annotated training data.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' pos_setParams(directory=datapath)
pos_setParams <- function(directory,
                          when_use_f0=10,
                          when_use_f0_first=10,
                          verbose=TRUE,
                          plot=TRUE) {

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
      print(paste0('Testing burst detection parameter settings for ', wav_list[snd]))
    }

    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

    i = 1

    for (ci in c(5, 10, 15, 20)) {

      for (rp in c(9, 12, 15, 18, 21)) {
        tmp <- positiveVOT(sig[,1], fs, vo_method='acf', closure_interval = ci,
                           release_param=rp,
                           burst_only=T, plot=F)

        pred$ci[i] <- ci
        pred$rp[i] <- rp
        pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$rel)
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
  good_ci <- pred$ci[winner]
  good_rp <- pred$rp[winner]

  finetune <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Finetuning burst detection parameter settings for ', wav_list[snd]))
    }
    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

    i = 1

    for (ci in c(good_ci-2, good_ci-1, good_ci, good_ci+1, good_ci+2)) {

      for (rp in c(good_rp-1, good_rp-0.5, good_rp, good_rp+0.5, good_rp+1)) {
        tmp <- positiveVOT(sig[,1], fs, vo_method='acf', closure_interval = ci,
                           release_param=rp,
                           burst_only=T, plot=F)

        finetune$ci[i] <- ci
        finetune$rp[i] <- rp
        finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$rel)
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
  n_bad <- length(which(finetune[winner,3:(ncol(finetune)-1)] >
                          (fs/1000)*when_use_f0_first))
  ci_winner <- finetune$ci[winner]
  rp_winner <- finetune$rp[winner]

  if (n_bad > (length(wav_list)/2)) {
    f0_first <- TRUE
    ci_winner <- NULL

    if (verbose) {
      print(paste0('Poor results achieved for most of the training material, ',
                   'retrying with f0_first=TRUE.'))
    }
    pred <- list()

    for (snd in 1:length(wav_list)){
      if (verbose) {
        print(paste0('Testing parameters with f0_first=TRUE for ', wav_list[snd]))
      }
      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i = 1

      for (rp in c(9, 12, 15, 18, 21)) {
        tmp <- positiveVOT(sig[,1], fs, release_param=rp,
                           burst_only=T, f0_first=T, plot=F)

        pred$rp[i] <- rp
        pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$rel)
        i <- i+1
      }

    }

    pred <- as.data.frame(pred)
    pred[is.na(pred)] <- 0.2*fs
    pred$diff_mean <- rowMeans(pred[2:ncol(pred)])

    if (verbose) {
      print(paste0('On average, the selected burst detection parameters after a ',
                   'first pass with f0_first=TRUE ',
                   'agree with the training data within a margin of ',
                   round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(pred$diff_mean)
    good_rp <- pred$rp[winner]

    finetune <- list()

    for (snd in 1:length(wav_list)){
      if (verbose) {
        print(paste0('Finetuning parameter settings with f0_first=TRUE for ',
                     wav_list[snd]))
      }
      wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
      tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

      sig <- wav[['sig']]
      fs <- wav[['fs']]
      tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

      i = 1

      for (rp in c(good_rp-1, good_rp-0.5, good_rp, good_rp+0.5, good_rp+1)) {
        tmp <- positiveVOT(sig[,1], fs,
                           release_param=rp,
                           burst_only=T, f0_first=T, plot=F)

        finetune$rp[i] <- rp
        finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$rel)
        i <- i+1
      }

    }

    finetune <- as.data.frame(finetune)
    finetune[is.na(finetune)] <- 0.2*fs
    finetune$diff_mean <- rowMeans(finetune[2:ncol(finetune)])

    if (verbose) {
      print(paste0('On average, the selected burst detection parameters after finetuning ',
                   'with method f0_first=TRUE ',
                   'agree with the training data within a margin of ',
                   round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(finetune$diff_mean)[1]
    rp_winner <- finetune$rp[winner]

  } else {
    f0_first <- FALSE
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
    tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

    i = 1

    for (gran in c(0.6, 0.8, 1, 1.2, 1.4)) {
      for (p in c(0.75, 0.8, 0.85, 0.9)) {
        tmp <- positiveVOT(sig[,1], fs, vo_method='acf', closure_interval=ci_winner,
                           release_param=rp_winner, f0_first=f0_first,
                           vo_granularity = gran, vo_param = p,
                           plot=F)

        pred$gran[i] <- gran
        pred$p[i] <- p
        pred[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$vo)
        i <- i+1
      }
    }
  }

  pred <- as.data.frame(pred)
  pred$diff_mean <- rowMeans(pred[3:ncol(pred)])

  if (verbose) {
    print(paste0('On average, the selected voicing onset parameters after a ',
                 'first pass ',
                 'agree with the training data within a margin of ',
                 round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
  }
  winner <- which.min(pred$diff_mean)
  good_gran <- pred$gran[winner]
  good_p <- pred$p[winner]

  finetune <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Finetuning voicing onset parameter settings for ', wav_list[snd]))
    }
    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

    i = 1

    for (gran in c(good_gran-0.1, good_gran-0.05, good_gran, good_gran+0.05,
                   good_gran+0.1)) {

      for (p in c(good_p-0.02, good_p-0.01, good_p, good_p+0.01, good_p+0.02)) {
        tmp <- positiveVOT(sig[,1], fs, vo_method='acf', closure_interval = ci_winner,
                           release_param=rp_winner, f0_first=f0_first,
                           vo_granularity = gran, vo_param = p,
                           plot=F)

        finetune$gran[i] <- gran
        finetune$p[i] <- p
        finetune[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$vo)
        i <- i+1
      }
    }
  }

  finetune <- as.data.frame(finetune)
  finetune$diff_mean <- rowMeans(finetune[3:ncol(finetune)])

  if (verbose) {
    print(paste0('On average, the selected voicing onset parameters after finetuning ',
                 'agree with the training data within a margin of ',
                 round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
  }
  winner <- which.min(finetune$diff_mean)[1]
  p_winner <- finetune$p[winner]
  gran_winner <- finetune$gran[winner]
  acf_results <- min(finetune$diff_mean) / (fs/1000)

  if (acf_results > when_use_f0) {
    if (verbose) {
      print(paste0('Average agreement with training data below ,',
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
      tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

      i = 1

      for (wl in c(20, 30, 40, 50, 60)) {
        for (acf in c(0.4, 0.5, 0.6, 0.7, 0.8)) {
          tmp <- positiveVOT(sig[,1], fs, vo_method='f0', closure_interval=ci_winner,
                             release_param=rp_winner, f0_first=f0_first,
                             f0_wl = wl, f0_minacf = acf,
                             plot=F)

          pred$wl[i] <- wl
          pred$acf[i] <- acf
          pred[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$vo)
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
      tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

      i = 1

      for (wl in c(good_wl-5, good_wl-2.5, good_wl, good_wl+2.5, good_wl+5)) {
        for (acf in c(good_acf-0.05, good_acf-0.025, good_acf, good_acf+0.025, good_acf+0.05)) {
          tmp <- positiveVOT(sig[,1], fs, vo_method='f0', closure_interval = ci_winner,
                             release_param=rp_winner, f0_first=f0_first,
                             f0_wl=wl, f0_minacf=acf,
                             plot=F)

          finetune$wl[i] <- wl
          finetune$acf[i] <- acf
          finetune[[wav_list[snd]]][i] <- abs(tg_t2 - tmp$vo)
          i <- i+1
        }
      }
    }

    finetune <- as.data.frame(finetune)
    finetune$diff_mean <- rowMeans(finetune[3:ncol(finetune)])

    if (verbose) {
      print(paste0('On average, the selected voicing parameters after finetuning ',
                   'with method=f0 agree with the training data within a margin of ',
                   round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
    }
    winner <- which.min(finetune$diff_mean)[1]
    wl_winner <- finetune$wl[winner]
    acf_winner <- finetune$acf[winner]
    f0_results <- min(finetune$diff_mean) / (fs/1000)

    if (f0_results < acf_results) {
      method <- 'f0'
    } else {
      method <- 'acf'
    }
  } else {
    method <- 'acf'
    wl_winner <- NULL
    acf_winner <- NULL
  }

  pl <- list(
    closure_interval = ci_winner,
    vo_granularity = gran_winner,
    vo_param = p_winner,
    release_param = rp_winner,
    f0_first = f0_first,
    vo_method = method,
    f0_wl = wl_winner,
    f0_minacf = acf_winner
  )

  if (plot) {
    plot_training_diff(directory, pl)
  }

  return(pl)

}
