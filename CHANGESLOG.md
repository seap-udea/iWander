
Changes log
-----------

The first version of the software was that published along the
submission of the paper to the arXiv. We call it the 0.1 version.

Since that initial version we have performed improvements on the
package either suggested by the reviewers of the paper, by colleagues
or by the community.  These are the most important
ones in the available versions of the manuscript.

### v0.1

This is the initial version.

### v0.2

Post-submission version.  These are the changes introduced to the
package:

* A new global configuration file ``iwander.conf`` was created that is
  common to all programs.  This configuration file now include the
  name of the analysed wanderer, eg. Oumuamua, Voyager1, etc.

* The output files from the main programs have as suffix the name of
  the wanderer, eg. wanderer-Oumuamua.csv, candidates-Oumuamua.csv.

* We have included new rules to makefile: all, cleandata, cleanall.

* Include elements of Voyager 1 and Voyager 2.

* New  scripts: 
  
  * progenitors.py: create a table of progenitor candidates in format
    markdown.

* The package allows now to compute the progenitor candidates and the
  future close encounters.  This is achieved by changing the sign of
  the ``duration```vairble in iwander.conf.

* List of candidates generated with progenitors.py, are different if
  they will be progenitor candidates or future close encounter
  candidates.

* ``probability.exe`` now calculates the 5%,50% and 95% percentiles of
  the critical quantities, tmin, dmin and vrel.

* ``probability.exe`` now uses the HIP2 properties in the database
  that are updated astrometric properties.

* We have added a new program to the package: ``reconstruct.cpp``.
  This program reconstruct the trajectories of the test particles and
  the progenitor candidates.

* The latest JPL solution for the Oumumua orbit in the Solar System
  (JPL-13 orbit) has been included.

* The solar galactocentric distance and its height above the mid
  galactic plane has been changed to updated values of RSUN = 8.2 kpc
  and ZSUN = 17 pc (Karim & Mamajek, 2016).

* The truncation radius has been reduced to 10^4 AU.

### v0.3

Version published after revision 1 of the paper:

* A new script intended to create the configuration file for a given
  object has been developed. The script is called JPL2iWander.py. A
  file called "<object>.jpl" must be provided and the script convert
  it to an "<object>.conf" file that is load in wanderer.conf.

* Configuration files for several objects has been added in the
  directory objects.  In the same directory the CANDIDATES tables will
  be stored.

