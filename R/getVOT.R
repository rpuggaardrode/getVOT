getVOT <- function(sound, sr, plot=TRUE,
                   closure_interval = 10,
                   release_param = 20,
                   vo_granularity = 1,
                   vo_param = 0.85) {

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
    sq_vo <- seq(from=rel+(step*5), to=rel+(step*200), by=(step*vo_granularity))
    mu_acf <- c()
    i <- 1
    for (s in sq_vo) {
      acf <- stats::acf(sound[s:(s+(step*vo_granularity)-1)], plot=F, na.action=na.pass)
      mu_acf[i] <- mean(acf$acf)
      i <- i+1
    }

    # hi_acf <- which(mu_acf > acf_threshold)
    hi_acf <- which(mu_acf > max(mu_acf, na.rm=T)*vo_param)
    if (is.na(hi_acf[2])) {
      vo <- (rel + (hi_acf[1]*(step*vo_granularity)) + step)
    } else {
      vo <- (rel + (hi_acf[2]*(step*vo_granularity)) + step)
    }

    vot <- round((vo-rel)/sr, 4) * 1000

    if(plot){
      plot(y=sound, x=seq(0, dur, length.out=length(sound)), type='l',
           xlab='Time (s)',
           ylab='Amplitude (dB)')
      abline(v=rel/sr, col='red')
      abline(v=vo/sr, col='red')
    }

    return(list(
      rel = rel,
      vo = vo,
      vot = vot
    ))
  }

  # f0 <- phonTools::pitchtrack(sound[rel:length(sound)], fs=sr, show=F,
  #                             windowlength=wl, minacf=minacf) #timestep?
  # vo <- round(f0$time[1]) * step + rel

  # sq_vo <- seq(from=rel+step, to=length(sound), by=step)
  # min_amp <- c()
  # i <- 1
  # for (s in sq_vo) {
  #   min_amp[i] <- min(sound[s:(s+step-1)])
  #   i <- i+1
  # }
  #
  # gulf_size <- min(sound)/2
  # gulf <- which(min_amp < gulf_size)
  # vo <- (rel + (gulf[1]*step))-step


}
