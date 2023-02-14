addVOT2emuDB <- function(emuDB, level, seg_list,
                         new_level_name = 'vot',
                         verbose = TRUE,
                         closure_interval = 10,
                         release_param = 15,
                         vo_method = 'acf',
                         vo_granularity = 1,
                         vo_param = 0.85,
                         f0_wl = 30, f0_minacf = 0.5){

  if (verbose) {
    print('Rewriting JSON files with new level definition')
  }
  emuR::add_levelDefinition(emuDB, new_level_name, type='SEGMENT', verbose=verbose)

  if ('character' %in% class(seg_list)) {
    if (verbose) {
      print('Finding segment intervals')
    }

    str <- paste0('ˈ', seg_list)
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

  #Looks like this bit doesn't actually do anything,
  #because when the mediafile_samples are imported,
  #emuR doesn't look at sample_start but looks at some of the other
  #vague identifiers in the tibble. Shucks.
  sl$sample_start <- sl$sample_start - (0.015*sl$sample_rate)
  sl$sample_start <- as.integer(sl$sample_start)

  if (verbose) {
    print('Importing WAV data')
  }

  suppressWarnings(
    all_wav <- get_trackdata(emuDB, sl, ssffTrackName='MEDIAFILE_SAMPLES',
                             verbose=verbose))

  if (verbose) {
    print('Searching for VOTs')
  }

  len <- length(unique(all_wav$sl_rowIdx))
  itc <- data.frame(session = rep(NA, len*2),
                    bundle = rep(NA, len*2),
                    level = rep('vot', len*2),
                    start_item_seq_idx = rep(NA, len*2),
                    attribute = rep('vot', len*2),
                    labels = rep(NA, len*2))

  for (i in 1:len) {
    one_wav <- all_wav[which(all_wav$sl_rowIdx == i),]
    vot <- getVOT(one_wav$T1, sl$sample_rate[i], plot=F,
                  closure_interval,
                  release_param, vo_method,
                  vo_granularity, vo_param, f0_wl, f0_minacf)
    t_start <- 15 + (sl$sample_start[i] / sl$sample_rate[i] * 1000) +
      (vot$rel / sl$sample_rate[i] * 1000)
    t_end <- 15 + (sl$sample_start[i] / sl$sample_rate[i] * 1000) +
      (vot$vo / sl$sample_rate[i] * 1000)

    itc$session[(i*2-1):(i*2)] <- sl$session[i]
    itc$bundle[(i*2-1):(i*2)] <- sl$bundle[i]
    itc$start_item_seq_idx[i*2-1] <- t_start
    itc$start_item_seq_idx[i*2] <- t_end
    itc$labels[i*2-1] <- 'vot'
    itc$labels[i*2] <- ''
  }

  if (verbose) {
    print('Rewriting JSON files with predicted VOT segmentation')
  }

  emuR::create_itemsInLevel(emuDB, itc, verbose=verbose)
  co <- emuR::get_levelCanvasesOrder(emuDB, 'default')
  co <- c(co, 'vot')
  emuR::set_levelCanvasesOrder(emuDB, 'default', co)

}
