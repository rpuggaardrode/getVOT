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
#' @param when_use_f0 Numeric; default is `10`. Per default, voicing onset is
#' predicted with [negativeVOT] using `method='acf'`. This is the preferred method
#' because it is significantly faster than the alternative `method='f0'`. For
#' this reason, parameters with `method='f0'` are only tested if the average
#' prediction error with `method='acf'` is above a certain threshold.
#' `when_use_f0` sets this threshold; i.e. by default, parameters for
#' `method='f0'` are tested when the average prediction error with `method='acf'`
#' is above 10 ms.
#' @param verbose Logical; should messages be printed in the console indicating
#' the progress? Default is `TRUE`.
#' @param plot Logical; should predicted VOT using the resulting parameters
#' be plotted against the hand-annotated VOT?
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
neg_setParams <- function(directory,
                          when_use_f0=10,
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
          pred[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
          i <- i+1
        }
      }
    }

  }

  pred <- as.data.frame(pred)
  pred$diff_mean <- rowMeans(pred[4:ncol(pred)])

  if (verbose) {
    print(paste0('On average, the selected parameters for voicing onset after a first pass ',
                 'agree with the training data within a margin of ',
                 round(min(pred$diff_mean) / (fs/1000), 3), ' ms'))
  }

  winner <- which.min(pred$diff_mean)
  good_ci <- pred$ci[winner]
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
    tg_t1 <- round(tg[['vot']][['t1']][2] * fs)

    i <- 1

    for (ci in c(good_ci-2, good_ci-1, good_ci, good_ci+1, good_ci+2)) {

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
          finetune[[wav_list[snd]]][i] <- abs(tg_t1 - tmp$vo)
          i <- i+1
        }
      }
    }

  }

  finetune <- as.data.frame(finetune)
  finetune$diff_mean <- rowMeans(finetune[4:ncol(finetune)])

  if (verbose) {
    print(paste0('On average, the selected voicing onset parameters after finetuning ',
                 'agree with the training data within a margin of ',
                 round(min(finetune$diff_mean) / (fs/1000), 3), ' ms'))
  }

  winner <- which.min(finetune$diff_mean)
  ci_winner <- finetune$ci[winner]
  gran_winner <- finetune$gran[winner]
  p_winner <- finetune$p[winner]
  acf_results <- min(finetune$diff_mean) / (fs/1000)

  if (acf_results > when_use_f0) {
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
    wl_winner <- 50
    acf_winner <- 0.5
  }

  pred <- list()

  for (snd in 1:length(wav_list)){
    if (verbose) {
      print(paste0('Testing accuracy of different burst detection methods ',
                   wav_list[snd]))
    }

    wav <- rPraat::snd.read(paste0(directory, '/', wav_list[snd]))
    tg <- rPraat::tg.read(paste0(directory, '/', tg_list[snd]))

    sig <- wav[['sig']]
    fs <- wav[['fs']]
    tg_t2 <- round(tg[['vot']][['t1']][3] * fs)

    for (rm in c('transient', 'amplitude')) {
      tmp <- negativeVOT(sig[,1], fs, vo_method=method,
                         f0_wl=wl_winner, f0_minacf=acf_winner,
                         closure_interval=ci_winner,
                         vo_granularity=gran_winner,
                         vo_param=p_winner,
                         rel_method=rm,
                         plot=F)

      pred[[rm]][snd] <- abs(tg_t2 - tmp$rel)
    }
  }

  if (mean(pred$transient) > mean(pred$amplitude)) {
    rm <- 'amplitude'
  } else {
    rm <- 'transient'
  }

  if (verbose) {
    print(paste0('On average, the select burst detection method agree',
                 'with the training data within a margin of ',
                 round(mean(pred[[rm]]) / (fs/1000), 3), ' ms'))
  }

  pl <- list(
    closure_interval = ci_winner,
    vo_granularity = gran_winner,
    vo_param = p_winner,
    vo_method = method,
    f0_wl = wl_winner,
    f0_minacf = acf_winner,
    rel_method=rm
  )

  if (plot) {
    plot_training_diff(directory, pl)
  }

  return(pl)

}
