VOT2newTG <- function(folder, verbose=TRUE,
                      closure_interval = 10,
                      release_param = 15,
                      vo_method = 'acf',
                      vo_granularity = 1,
                      vo_param = 0.85,
                      f0_wl = 30, f0_minacf = 0.5) {

  fls <- list.files(folder)
  dir.create(paste0(folder, '/tg'))

  for (f in fls){
    flname <- unlist(strsplit(f, '[.]'))[1]
    fltype <- unlist(strsplit(f, '[.]'))[2]

    if (fltype != 'wav') {
      if (verbose) {
        print(paste0(f, ' was not processed as this function only works with
                   sound files with the .wav extension'))
      }
    } else {
      snd <- rPraat::snd.read(f)
      vot <- getVOT(snd$sig, snd$fs, plot=F, closure_interval,
                    release_param, vo_method,
                    vo_granularity, vo_param, f0_wl, f0_minacf)

      tg <- rPraat::tg.createNewTextGrid(0, snd$duration)
      tg <- rPraat::tg.insertNewIntervalTier(tg, newTierName='vot')
      tg <- rPraat::tg.insertInterval(tg, 'vot',
                              tStart=vot$rel/snd$fs,
                              tEnd=vot$vo/snd$fs,
                              label=as.character(vot$vot))

      rPraat::tg.write(tg, paste0(folder, '/tg/', flname, '.TextGrid'))

      if (verbose) {
        print(paste0(f, ' processed'))
      }
    }
  }
}
