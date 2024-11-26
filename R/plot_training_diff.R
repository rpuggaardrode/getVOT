#' Plot hand-annotated voice onset time and predicted voice onset time together
#'
#' Used to plot hand-annotated voice onset time ("training data") together with
#' voice onset time that has been predicted using optimized parameters based on
#' that training data. This is practical in order to check the performance of
#' the selected parameters. Predicted voice onset time is shown in red, and
#' hand-annotated voice onset time is shown in blue.
#'
#' @param directory A string giving the location of a directory with WAV/TextGrid
#' pairs, previously used to optimize parameters with [neg_setParams] or
#' [pos_setParams].
#' @param params_list A named list of optimized parameters returned by
#' [neg_setParams] or [pos_setParams].
#' @param cat Logical; default is `FALSE`. If there are empty frames in the
#' resulting plot, should they be filled with pictures of my very nice cat
#' Kipawsky?
#'
#' @return Called for side effects; produces a plot.
#' @export
#'
#' @examples
#' datapath <- system.file('extdata/vl', package='getVOT')
#' p <- neg_setParams(directory=datapath, plot=FALSE)
#' plot_training_diff(directory=datapath, params_list=p)
plot_training_diff <- function(directory,
                               params_list,
                               cat=FALSE) {

  p <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(p))

  wavs <- list.files(directory, '*.wav')
  tgs <- list.files(directory, '*.TextGrid')

  num <- length(wavs) + 1
  frames <- ceiling(sqrt(num))
  if (num > frames*(frames-1)) {
    plot_dim <- c(frames, frames)
  } else {
    plot_dim <- c(frames, (frames-1))
  }

  graphics::par(mfrow=plot_dim, mai=c(0.55, 0.55, 0.15, 0.15))

  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  graphics::legend(x = 'center', y = 'center',
                   legend = c('Predicted', 'Hand-annotated'),
                   col = c('red', 'blue'), lty = c('solid', 'dashed'),
                   inset = 0, cex = 1.2, lwd = 2)

  for (s in 1:(num-1)) {
    snd <- rPraat::snd.read(paste0(directory, '/', wavs[s]))
    tg <- rPraat::tg.read(paste0(directory, '/', tgs[s]))
    sig <- snd[['sig']]
    fs <- snd[['fs']]

    if ('f0_first' %in% names(params_list)) {
      tmp <- positiveVOT(sig[,1], fs, params_list=params_list)
      if (tmp$vot == "NA") {
        plot(y=sig, x=seq(0, length(sig)/fs, length.out=length(sig)), type='l',
             xlab='Time (s)', ylab='Amplitude')
      }
    } else {
      tmp <- negativeVOT(sig[,1], fs, params_list=params_list)
      if (tmp$vot == "NA") {
        plot(y=sig, x=seq(0, length(sig)/fs, length.out=length(sig)), type='l',
             xlab='Time (s)', ylab='Amplitude')
      }
    }

    graphics::abline(v=tg[['vot']][['t1']][2], col='blue', lty='dashed', lwd=2)
    graphics::abline(v=tg[['vot']][['t1']][3], col='blue', lty='dashed', lwd=2)
  }


  if (cat) {
    spots <- plot_dim[1]*plot_dim[2]
    if (num < spots) {
      cat_pic <- jpeg::readJPEG(system.file('extdata/cats/cropcat.jpg', package='getVOT'))
      plot(1, type='n', ylab='cat', xlab='cat', xaxt='n', yaxt='n')
      graphics::rasterImage(cat_pic, 0.6, 0.6, 1.4, 1.4)
    }
    if (num < (spots-1)) {
      cat_pic <- jpeg::readJPEG(system.file('extdata/cats/cropcat2.jpg', package='getVOT'))
      plot(1, type='n', ylab='cat', xlab='cat', xaxt='n', yaxt='n')
      graphics::rasterImage(cat_pic, 0.6, 0.6, 1.4, 1.4)
    }
  }

  graphics::par(p)
}
