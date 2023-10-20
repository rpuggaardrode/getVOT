#' Predict voice onset time and add results to EMU database
#'
#' Predicts voice onset time for all intervals in an EMU database that
#' match a pattern provided by the user. Adds an additional level definition
#' to all annotations in the EMU database giving the location of predicted
#' voice onset time (i.e. pre-voicing or stop release).
#'
#' @param emuDB The handle of an EMU database which is already loaded in R.
#' @param seg_list A segment list returned by [emuR::query] containing stops
#' of interest. Alternatively one or more strings giving the characters used to
#' annotate stops in `level`. In this case, does not have to be a perfect match;
#' the function will
#' search for all intervals where the first character matches one of the strings
#' in `seg_list`. If any intervals contain the IPA stress symbol, this symbol
#' will be ignored. [emuR::query] is then called under the hood.
#' @param level A string giving the name of the annotation level in `emuDB`
#' to search for instances of segments given in `seg_list`. Default is `NULL`;
#' not used if a `seg_list` returned by [emuR::query] is provided.
#' @param new_level_name A string given the name of the new annotation level
#' used to store predicted voice onset time. Default is `vot`.
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
#' @seealso This function is intended for users of the EMU-SDMS database system
#' and the [emuR] package. Sister functions which are based on the acoustic
#' processing software Praat are also available in [addVOT2TG] and [VOT2newTG],
#' and these are significantly faster and less computationally expensive; if
#' possible, it is recommended to use these instead.
#' Voice onset time is predicted using the [getVOT] function, which
#' calls either the [positiveVOT] or [negativeVOT] function (or both) for a
#' given token.
#' @export
#'
#' @examples
#' #later
addVOT2emuDB <- function(emuDB, seg_list, level=NULL,
                         new_level_name = 'vot',
                         verbose = TRUE,
                         sign=c('positive', 'negative'),
                         pos_params_list=NULL,
                         neg_params_list=NULL){

  if (verbose) print('Rewriting JSON files with new level definition')
  emuR::add_levelDefinition(emuDB, new_level_name, type='SEGMENT', verbose=verbose)

  if ('character' %in% class(seg_list)) {
    if (verbose) print('Finding segment intervals')
    str <- paste0('\\u02c8', seg_list)
    labs <- c(seg_list, str)

    sl <- emuR::query(emuDB, paste0(level, " =~ '^", labs[1], ".*'"))
    for (lab in labs[-1]) {
      sl <- suppressMessages(
        dplyr::full_join(sl,
                         emuR::query(emuDB, paste0(level, " =~ '^", lab, ".*'"))))
    }
  } else {
    sl <- seg_list
  }

  len <- nrow(sl)

  itc <- data.frame(session = rep(NA, len*2),
                    bundle = rep(NA, len*2),
                    level = rep(new_level_name, len*2),
                    start_item_seq_idx = rep(NA, len*2),
                    attribute = rep(new_level_name, len*2),
                    labels = rep(NA, len*2))

  dbloc <- emuDB$basePath

  if (verbose) print('Predicting VOT')

  for (i in 1:len) {
    fn <- paste0(dbloc, '/', sl$session[i], '_ses/', sl$bundle[i], '_bndl/',
                 sl$bundle[i], '.wav')
    if (sl$sample_start[i] - (0.015*sl$sample_rate[i]) > 0) {
      t1 <- round(sl$sample_start[i] - (0.015*sl$sample_rate[i]))
    } else {
      t1 <- sl$sample_start[i]
    }

    snd <- rPraat::snd.read(fn, from=t1, to=sl$sample_end[i])
    vot <- getVOT(snd$sig[,1], sl$sample_rate[i], sign=sign,
                  neg_params_list, pos_params_list)

    if (vot$vot > 0) {
      t_start <- (t1 / sl$sample_rate[i] * 1000) +
        (vot$rel / sl$sample_rate[i] * 1000)
      t_end <- (t1 / sl$sample_rate[i] * 1000) +
        (vot$vo / sl$sample_rate[i] * 1000)
    } else {
      t_start <- (t1 / sl$sample_rate[i] * 1000) +
        (vot$vo / sl$sample_rate[i] * 1000)
      t_end <- (t1 / sl$sample_rate[i] * 1000) +
        (vot$rel / sl$sample_rate[i] * 1000)
    }

    itc$session[(i*2-1):(i*2)] <- sl$session[i]
    itc$bundle[(i*2-1):(i*2)] <- sl$bundle[i]
    itc$start_item_seq_idx[i*2-1] <- t_start
    itc$start_item_seq_idx[i*2] <- t_end
    itc$labels[i*2-1] <- 'vot'
    itc$labels[i*2] <- ''
  }

  if (verbose) print('Rewriting JSON files with predicted VOT segmentation')

  emuR::create_itemsInLevel(emuDB, itc, verbose=verbose)
  co <- emuR::get_levelCanvasesOrder(emuDB, 'default')
  co <- c(co, new_level_name)
  emuR::set_levelCanvasesOrder(emuDB, 'default', co)

}
