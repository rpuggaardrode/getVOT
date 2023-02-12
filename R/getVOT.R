getVOT <- function(sound, sr, plot=TRUE,
                   closure_interval = 10,
                   release_param = 15,
                   vo_method = 'acf',
                   vo_granularity = 1,
                   vo_param = 0.85,
                   f0_wl = 30, f0_minacf = 0.5) {

  ci <- closure_interval/1000

  dur <- length(sound)/sr
  sqlen <- ceiling(length(sound)/2 / sr / ci)
  sq_clo <- seq(from=sr*ci, to=sqlen*sr*ci, by=sr*ci)

  mean_amp <- c()
  i <- 1
  for (s in sq_clo) {
    mean_amp[i] <- mean(abs(sound[(s-(sr*ci)+1):s]))
    i <- i+1
  }
  clo <- (which(mean_amp==min(mean_amp)) - 0.5) * sr * ci

  step <- 0.001*sr
  sq_rel <- seq(from=clo[1], to=length(sound), by=step)
  max_amp <- c()
  i <- 1
  for (s in sq_rel) {
    max_amp[i] <- max(sound[s:(s+step-1)])
    i <- i+1
  }

  spike_size <- max(sound)/release_param
  spike <- which(max_amp > spike_size)
  rel <- (clo[1] + ((spike[1])*step))-step

  if (is.na(rel)) {
    return(list(
      rel = clo,
      vo = clo+(0.02*sr),
      vot = 'NA'
    ))
  } else {
    if (!(vo_method %in% c('acf', 'f0'))) {
      stop('Legal parameters for vo_method are acf and f0.')
    }

    if (vo_method == 'acf') {
      sq_vo <- seq(from=rel+(step*5), to=rel+(step*200), by=(step*vo_granularity))
      mu_acf <- c()
      i <- 1
      for (s in sq_vo) {
        acf <- stats::acf(sound[s:(s+(step*vo_granularity)-1)], plot=F, na.action=na.pass)
        mu_acf[i] <- mean(acf$acf)
        i <- i+1
      }

      hi_acf <- which(mu_acf > max(mu_acf, na.rm=T)*vo_param)
      if (is.na(hi_acf[2])) {
        vo <- (rel + (hi_acf[1]*(step*vo_granularity)) + step)
      } else {
        vo <- (rel + (hi_acf[2]*(step*vo_granularity)) + step)
      }


    } else if (vo_method == 'f0') {
      f0 <- phonTools::pitchtrack(sound[rel:length(sound)], fs=sr, show=F,
                                  windowlength=f0_wl, minacf=f0_minacf)
      vo <- round(f0$time[1]) * step + rel

      if (is.na(vo)) {
        return(list(
          rel = rel,
          vo = rel+(0.02*sr),
          vot = 'NA'
        ))
      }
    }

    vot <- round((vo-rel)/sr, 4) * 1000

    if(plot){
      plot(y=sound, x=seq(0, dur, length.out=length(sound)), type='l',
           xlab='Time (s)',
           ylab='Amplitude')
      abline(v=rel/sr, col='red')
      abline(v=vo/sr, col='red')
    }

    return(list(
      rel = rel,
      vo = vo,
      vot = vot
    ))
  }

}

par(mfrow=c(3,3))

for (f in fls){
  snd <- rPraat::snd.read(f)
  getVOT(snd$sig, snd$fs)
}
