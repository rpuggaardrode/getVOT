
<!-- README.md is generated from README.Rmd. Please edit that file -->

# getVOT

`getVOT` is a developmental R package with functions that aim to predict
positive or negative voice onset time from bare or annotated sound
files. VOT prediction itself is done by loading (partial) sound files
into R and using largely base R functions to predict the usual VOT
landmarks from a vector of audio samples. Options are available to
generate [Praat TextGrids](https://fon.hum.uva.nl/praat) or add a tier
to existing TextGrids with predicted VOT, or to add an annotation level
with predicted VOT to an existing [EMU
database](https://github.com/IPS-LMU/emuR) (although the latter
functionality is significantly slower).

The goal is to partially make automatic, or at least aid in, the job of
annotating VOT. There are other and more sophisticated algorithms
available for predicting VOT, including
[AutoVOT](https://github.com/mlml/autovot/tree/master) and
[Dr.VOT](https://github.com/MLSpeech/Dr.VOT). The job of `getVOT` is not
to replace or even necessarily outperform these. AutoVOT only predicts
positive VOT and scales best with a fair amount of training data; Dr.VOT
predicts both positive and negative VOT, but can be a bit daunting to
set up and difficult to manipulate; both can only be used on Unix-style
operating systems. `getVOT` should be easy to set up and get started
with, should be compatible with all operating systems, and should scale
fairly well with little to no training data.\* This makes `getVOT`
optimal for smaller data sets where providing suitable training data is
difficult.

A big motivation for creating `getVOT` was also as a coding exercise. I
felt that it *should* be possible to find the regular VOT landmarks
using a relatively simple algorithm with few bells and whistles, and I
have personally been rather happy with the results. Because I know
exactly what `getVOT` does under the hood, it’s also usually
straightforward to figure out why it fails when it fails.

\*This is a beta version of `getVOT`! If you find that any of these
claims do not match your experience, I’d be very happy to hear about it!
You can reach out at r.puggaard at phonetik.uni-muenchen.de.

## Installation

You can install the development version of `getVOT` from GitHub with:

``` r
#install.packages('devtools')
devtools::install_github('rpuggaardrode/getVOT')
```

## Workflow

I give a few examples below of how the main functions of `getVOT` work.
There is a longer tutorial on
[LingMethodsHub](https://lingmethodshub.github.io/content/R/getVOT-tutorial/),
which also gets into the nitty-gritty of how it works under the hood.

### No training data

Probably the two most important `getVOT` functions in daily use will be
`VOT2newTG` and `addVOT2TG`. As the names suggest, `VOT2newTG` will
generate a new TextGrid with predicted VOT where no annotations exist
already. In this case, sound files should be short and consist of only a
single word beginning with a stop. `addVOT2TG` adds a new tier to an
existing TextGrid with predicted VOT. In this case sound files can be
short, but the existing TextGrid should indicate the rough location of
stops.

`VOT2newTG` works like this:

``` r
VOT2newTG(directory='my_directory', sign='positive')
```

`my_directory` should be a directory of sound files. The `sign` argument
specifies whether `positive` or `negative` VOT should be predicted. This
argument is optional; the default is to use a simple algorithm to try to
predict whether VOT is positive or negative, but this algorithm is not
all that precise (yet).

The function will create a new subdirectory of `my_directory` containing
TextGrids with VOT annotations.

`addVOT2TG` works like this:

``` r
addVOT2TG(directory='my_directory', tg_tier='my_word_tier',
          seg_list=c('p', 't', 'k'), sign='positive')
```

In this case, `my_directory` should be a directory of sound-TextGrid
pairs. `tg_tier` gives the name of the TextGrid tier where the rough
locations of stops are indicated. `seg_list` is a list of stop symbols
to look for in `my_word_tier`; it does not need to be a precise match,
in this case VOT will be predicted for all intervals of `my_word_tier`
where the first symbol is `p`, `t`, or `k`, ignoring stress symbols.

### Training data

`VOT2newTG` and `addVOT2TG` also take the optional arguments
`pos_params_list` and `neg_params_list`; these are lists of some (rather
obscure) parameters that are used to tweak VOT prediction. Under the
hood, these parameters are passed on to the functions `positiveVOT` and
`negativeVOT`. These are the default arguments:

``` r
args(getVOT::positiveVOT)
#> function (sound, sr, closure_interval = 10, release_param = 15, 
#>     vo_method = "acf", vo_granularity = 1, vo_param = 0.85, f0_wl = 30, 
#>     f0_minacf = 0.5, burst_only = FALSE, f0_first = FALSE, plot = TRUE, 
#>     params_list = NULL) 
#> NULL
args(getVOT::negativeVOT)
#> function (sound, sr, vo_method = "acf", closure_interval = 10, 
#>     vo_granularity = 1.2, vo_param = 0.9, f0_wl = 50, f0_minacf = 0.5, 
#>     vo_only = FALSE, rel_method = "transient", plot = TRUE, params_list = NULL) 
#> NULL
```

More info about each of these can be found in the respective help files,
`?positiveVOT` and `?negativeVOT`. Mileage may vary with the defaults –
they were chosen because I have found that they usually work fairly
well, but no settings of these parameters will give equally good results
for all data. Parameter lists can be edited by hand, but `getVOT` also
offers the functions `pos_setParams` and `neg_setParams` which will help
optimize parameters for your data.

In order to use `pos_setParams` and `neg_setParams` you should provide a
few (say, 5-10) representative examples of sound file-TextGrid pairs,
where the TextGrid has a single tier with manual VOT annotation. The
`setParams` functions will then try out a range of different parameter
settings, and return a list with the settings that are found to minimize
the difference between predicted VOT and manually labeled VOT.

You simply need to pass the functions a directory, like so:

``` r
optimized_params <- pos_setParams(directory='my_training_data')
```

The resulting list `optimized_params` can then be passed onto functions
like `VOT2newTG` like so:

``` r
VOT2newTG(directory='my_directory', sign='positive',
          pos_params_list=optimized_params)
```
