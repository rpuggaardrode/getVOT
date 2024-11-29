# getVOT 0.2.1

* Made some slight improvements to the functionality for predicting whether
a stop is pre-voiced or not.

# getVOT 0.2.0

* The code has been extensively restructured and is now much more modular,
with most of the closure, release, and voicing onset detection methods stored
in separate functions.
* A new voicing onset detection method has been added using strength of 
harmonic excitation (`vo_method='soe'`). This is now the default for predicting
negative VOT, as it works quite well in the tests I've run and is much 
faster than all other methods. For this reason, the `when_avoid_soe` argument
has been added to `neg_setParams()`, specifying an error value where other
methods should be tried (default is 10 ms), since `neg_setParams()` can 
otherwise be very slow.
* A new burst detection method has been added, where the speech signal is
differenced before searching for sudden amplitude increases in the *differenced*
signal. This is now the default when looking for positive VOT, as it works
quite well in the tests I've run. 
* In order to speed up `negativeVOT` I've removed part of the functionality
for estimating whether a stop has prevoicing or not (the existing functionality
didn't work very well, and it wasn't really worth the speed costs). You can 
still specify that `getVOT()` should look for both positive and negative VOT,
but it'll likely perform even worse than before. I haven't found a solution to
this yet.
* The line widths in plots produced by `getVOT` has been increased, and a 
legend is now included when using `plot_training_diff()`.
