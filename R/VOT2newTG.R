#' Predict voice onset time and create TextGrid with the results
#'
#' Loops through a directory with WAV files, each containing a syllable with
#' an initial stop consonant. Predicts voice onset time and generates a new
#' Praat TextGrid file giving the location of the predicted voice onset time
#' (i.e. pre-voicing or stop release).
#'
#' @param directory A string giving the location of a directory with WAV files.
#' The results will be stored in a new subdirectory called `tg`.
#' @param verbose Logical; should messages be printed in the console indicating
#' the progress? Default is `TRUE`.
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
#' given token. See also the sister function [addVOT2TG] which adds
#' predicted voice onset time to an existing TextGrid where the rough location
#' of a stop has been annotated.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' VOT2newTG(datapath)
VOT2newTG <- function(directory, verbose=TRUE,
                      sign=c('positive', 'negative'),
                      neg_params_list=NULL,
                      pos_params_list=NULL) {

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

  fls <- list.files(directory)
  dir.create(paste0(directory, '/tg'))

  for (f in fls){
    flname <- unlist(strsplit(f, '[.]'))[1]
    fltype <- unlist(strsplit(f, '[.]'))[2]

    if (fltype != 'wav') {
      if (verbose) {
        print(paste0(f, ' was not processed as this function only works with
                   sound files with the .wav extension'))
      }
    } else {
      snd <- rPraat::snd.read(paste0(directory, '/', f))

      vot <- getVOT(snd$sig[,1], snd$fs, sign=sign, neg_params_list, pos_params_list)

      tg <- rPraat::tg.createNewTextGrid(0, snd$duration)
      tg <- rPraat::tg.insertNewIntervalTier(tg, newTierName='vot')
      if (vot$vot > 0) {
        tg <- rPraat::tg.insertInterval(tg, 'vot',
                                        tStart=vot$rel/snd$fs,
                                        tEnd=vot$vo/snd$fs,
                                        label='pos')
      } else {
        tg <- rPraat::tg.insertInterval(tg, 'vot',
                                        tStart=vot$vo/snd$fs,
                                        tEnd=vot$rel/snd$fs,
                                        label='neg')
      }

      rPraat::tg.write(tg,
                       paste0(directory, '/tg/', flname, '.TextGrid'),
                       format='text')

      if (verbose) {
        print(paste0(f, ' processed'))
      }
    }
  }
}
