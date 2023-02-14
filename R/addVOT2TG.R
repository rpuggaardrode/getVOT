addVOT2TG <- function(folder, tg_tier, seg_list,
                      verbose = TRUE,
                      new_tier_name = 'vot',
                      closure_interval = 10,
                      release_param = 15,
                      vo_method = 'acf',
                      vo_granularity = 1,
                      vo_param = 0.85,
                      f0_wl = 30, f0_minacf = 0.5){

  fls <- list.files(folder)
  dir.create(paste0(folder, '/vot'))

  wav <- fls[which(stringr::str_sub(fls, -3) == 'wav')]
  tg <- fls[which(stringr::str_sub(fls, -8) == 'TextGrid')]

  wav_flname <- strsplit(wav, '[.]')
  wav_flname <- unlist(purrr::map(wav_flname, 1))
  tg_flname <- strsplit(tg, '[.]')
  tg_flname <- unlist(purrr::map(tg_flname, 1))

  if (verbose) {
    if (any(wav_flname != tg_flname)) {
      print('There is no one-to-one correspondence between .wav files and
            .TextGrid files in this folder. Only .wav files with a corresponding
            .TextGrid file will be processed.')
    }
  }

  tg_flname <- tg_flname[which(tg_flname %in% wav_flname)]

  for (tgfl in tg_flname){

    tg <- rPraat::tg.read(paste0(tgfl, '.TextGrid'))
    tg <- rPraat::tg.duplicateTier(tg, tg_tier, newTierName='copy')
    tg$copy$label <- stringr::str_replace_all(tg$copy$label, 'ˈ', '')
    tg <- rPraat::tg.insertNewIntervalTier(tg, newTierName=new_tier_name)
    # sl_ints <- which(any(startsWith(tg$copy$label, seg_list))) #stuck at this ridiculous problem,
    # so currently everything in seg_list has to be one char.
    sl_ints <- which(stringr::str_sub(tg$copy$label, 0, 1) %in% seg_list)

    for (int in sl_ints) {
      if (tg$copy$t1[int] - 0.015 > 0) {
        t1 <- tg$copy$t1[int] - 0.015
      } else {
        t1 <- tg$copy$t1[int]
      }

      snd <- rPraat::snd.read(paste0(folder, '/', tgfl, '.wav'),
                              from=t1,
                              to=tg$copy$t2[int], units='seconds')
      vot <- getVOT(snd$sig, snd$fs, plot=F, closure_interval,
                    release_param, vo_method,
                    vo_granularity, vo_param, f0_wl, f0_minacf)
      tg <- rPraat::tg.insertInterval(tg, new_tier_name,
                                      tStart=t1 + vot$rel/snd$fs,
                                      tEnd=t1 + vot$vo/snd$fs,
                                      label='vot')
    }

    tg <- rPraat::tg.removeTier(tg, 'copy')
    rPraat::tg.write(tg, paste0(folder, '/vot/', tgfl, '.TextGrid'))
    if (verbose){
      print(paste0(tgfl, ' processed'))
    }

  }
}
