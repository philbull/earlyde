# earlyde
Fit early-ish dark energy models to background distance/expansion rate data.

Requirements
============
 * numpy, scipy, matplotlib
 * pyccl (for some test functions)
 * emcee


Notes
=====
 - For all of the results generated before 2020-06-15, there was a minus sign 
   bug in the definition of the Tracker model 'deltaw' parameter (if used). It 
   works out so that 'deltaw' in these runs is actually '-deltaw'. This bug has 
   now been fixed.
