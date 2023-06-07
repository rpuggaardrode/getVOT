#' Predict voice onset time
#'
#' Voice onset time and the location of the associated acoustic landmarks are
#' predicted for a sound snippet including a stop consonant.
#'
#' @param sound Numeric vector corresponding to a sound wave.
#' @param sr Integer; sample rate of `sound`.
#' @param sign A string giving the type of voice onset time to predict;
#' legal values are `positive`, `negative`, or both. Default is both, in which
#' case the function tries to predict whether voice onset time is positive or
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
#' @return A named list with the following elements:
#' * `vo` Gives the sample number where the onset of voicing is predicted.
#' * `rel` Gives the sample number where the stop burst is predicted.
#' * `vot` Gives the predicted voice onset time in ms.
#' @seealso This function calls either [positiveVOT], [negativeVOT] or both.
#' It is used under the hood for adding predicted VOT to a Praat TextGrid with
#' [VOT2newTG] or [addVOT2TG] or to an EMU database with [addVOT2emuDB].
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' snd <- rPraat::snd.read(paste0(datapath, '/1.wav'))
#' sig <- snd$sig[,1]
#' sr <- snd$fs
#'
#' getVOT(sound=sig, sr=sr)
getVOT <- function(sound, sr,
                   sign=c('positive', 'negative'),
                   neg_params_list=NULL,
                   pos_params_list=NULL) {

  if ('positive' %in% sign & 'negative' %in% sign) {
    pos_test <- positiveVOT(sound, sr, plot=F,
                            params_list=pos_params_list)

    if (pos_test$spike > 0.025 | pos_test$vot < 10) {
      neg_test <- negativeVOT(sound, sr, plot=F,
                              params_list=neg_params_list)

      if (neg_test$vo < pos_test$rel && -neg_test$vot < neg_test$voi_int/2) {
        return(neg_test)
      } else {
        return(pos_test[-4])
      }
    } else {
      return(pos_test[-4])
    }
  } else if ('positive' %in% sign) {
    tmp <- positiveVOT(sound, sr, plot=F, params_list=pos_params_list)
    return(tmp[-4])
  } else if ('negative' %in% sign) {
    tmp <- negativeVOT(sound, sr, plot=F, params_list=neg_params_list)
  }
}
