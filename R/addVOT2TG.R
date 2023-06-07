#' Predict voice onset time and add results to existing TextGrid
#'
#' Loops through a directory with pairs of WAV files and Praat TextGrids,
#' and predicts voice onset time for all intervals that match a pattern
#' provided by the user. Copies of the TextGrid files are generated with an
#' additional tier giving the location of the predicted voice onset time
#' (i.e. pre-voicing or stop release).
#'
#' @param directory A string giving the location of a directory with pairs of
#' WAV files and TextGrids.
#' The results will be stored in a new subdirectory called `tg`.
#' @param tg_tier String giving the name of a tier in the TextGrids containing
#' intervals with stops. Alternatively an integer giving the index of a tier
#' containing intervals with stops stops.
#' @param seg_list One or more strings giving the characters used to annotate
#' stops in `tg_tier`. Does not have to be a perfect match; the function will
#' search for all intervals where the first character matches one of the strings
#' in `seg_list`. If any intervals contain the IPA stress symbol, this symbol
#' will be ignored.
#' @param verbose Logical; should messages be printed in the console indicating
#' the progress? Default is `TRUE`.
#' @param new_tier_name String giving the name of the new tier containing
#' VOT predictions. Default is `vot`.
#' @param sign A string giving the type of voice onset time to predict;
#' legal values are `positive`, `negative`, or both. Default is both, in which
#' case [getVOT] tries to predict whether voice onset time is positive or
#' negative. The development of this function is still in the early stages,
#' so it is recommended to specify whether VOT is expected to be positive or
#' negative when possible.
#' @param neg_params_list Named list of parameters used to predict negative
#' voice onset time; see [negativeVOT] for more information. Default is `NULL`,
#' in which case default parameters are used, but note that the default
#' parameters do not necessarily scale particularly well.
#' @param pos_params_list Named list of parameters used to predict positive
#' voice onset time; see [positiveVOT] for more information. Default is `NULL`,
#' in which case default parameters are used, but note that the default
#' parameters do not necessarily scale particularly well.
#'
#' @seealso WAV files are imported and TextGrids are exported using the `rPraat`
#' package. Voice onset time is predicted using the [getVOT] function, which
#' calls either the [positiveVOT] or [negativeVOT] function (or both) for a
#' given token. See also the sister function [VOT2newTG] which generates new
#' TextGrids with predicted voice onset time for sound files that consist of
#' a single stop-initial syllable or word.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/add2tg', package='getVOT')
#' addVOT2TG(directory=datapath, tg_tier='KAN-MAU', seg_list='b')
addVOT2TG <- function(directory, tg_tier, seg_list,
                      verbose = TRUE,
                      new_tier_name = 'vot',
                      sign = c('positive', 'negative'),
                      pos_params_list=NULL,
                      neg_params_list=NULL
                      ){

  if (length(sign) == 1) {
    if (!(sign %in% c('positive', 'negative'))) {
      stop('The value of sign should be either positive, negative, or both')
    }
  } else if (length(sign) == 2) {
    if(!any(sign %in% c('positive', 'negative'))) {
      stop('The value of sign should be either positive, negative, or both')
    }
  } else {
    stop('The value of sign should be either positive, negative, or both')
  }

  sl_regex <- paste(paste0('^', seg_list), collapse='|')

  dir.create(paste0(directory, '/vot'))

  wav <- list.files(directory, '*.wav')
  tg <- list.files(directory, '*.TextGrid')

  wav_flname <- strsplit(wav, '[.]')
  wav_flname <- unlist(purrr::map(wav_flname, 1))
  tg_flname <- strsplit(tg, '[.]')
  tg_flname <- unlist(purrr::map(tg_flname, 1))

  if (verbose) {
    if (any(wav_flname != tg_flname)) {
      print(paste0('There is no one-to-one correspondence between .wav files and',
            '.TextGrid files in this directory. Only .wav files with a corresponding',
            '.TextGrid file will be processed.'))
    }
  }

  tg_flname <- tg_flname[which(tg_flname %in% wav_flname)]

  for (tgfl in tg_flname){

    tg <- rPraat::tg.read(paste0(directory, '/', tgfl, '.TextGrid'))
    tg <- rPraat::tg.duplicateTier(tg, tg_tier, newTierName='copy')
    tg$copy$label <- stringr::str_replace_all(tg$copy$label, '\\u02c8', '')
    tg <- rPraat::tg.insertNewIntervalTier(tg, newTierName=new_tier_name)
    sl_ints <- which(stringr::str_detect(tg$copy$label, sl_regex))

    for (int in sl_ints) {
      if (tg$copy$t1[int] - 0.015 > 0) {
        t1 <- tg$copy$t1[int] - 0.015
      } else {
        t1 <- tg$copy$t1[int]
      }

      snd <- rPraat::snd.read(paste0(directory, '/', tgfl, '.wav'),
                              from=t1,
                              to=tg$copy$t2[int], units='seconds')
      vot <- getVOT(snd$sig[,1], snd$fs, sign=sign, neg_params_list, pos_params_list)
      if (vot$vot > 0) {
        tg <- rPraat::tg.insertInterval(tg, new_tier_name,
                                        tStart=t1 + vot$rel/snd$fs,
                                        tEnd=t1 + vot$vo/snd$fs,
                                        label='pos')
      } else {
        tg <- rPraat::tg.insertInterval(tg, new_tier_name,
                                        tStart=t1 + vot$vo/snd$fs,
                                        tEnd=t1 + vot$rel/snd$fs,
                                        label='neg')
      }
    }

    tg <- rPraat::tg.removeTier(tg, 'copy')
    rPraat::tg.write(tg,
                     paste0(directory, '/vot/', tgfl, '.TextGrid'),
                     format='text')
    if (verbose){
      print(paste0(tgfl, ' processed'))
    }

  }
}
