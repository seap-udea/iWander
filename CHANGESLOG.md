
Changes log
-----------

The first version of the software was that published along the
submission of the paper to the arXiv. We call it the 0.1 version.

Since that initial version we have performed improvements on the
package either suggested by the reviewers of the paper, by colleagues
or by the community.  These are the most important
ones in the available versions of the manuscript.

### v0.1


Initial version.

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
