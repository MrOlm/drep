# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [1.2.0] - 2017-19-11
### fixed
- Scoring with regards to strain heterogeneity modified

## [1.1.3] - 2017-09-11
### fixed
- N50 calculation is now correct

## [1.1.2] - 2017-08-25
### Changed
- fixed proper pip sklean-learn (thanks Ben Woodcroft)
- added the blank folder test_backend/ so that the test suite will work

## [1.1.1] - 2017-07-26
### Changed
- added links to ISME publication in readme and documentation

## [1.1.0] - 2017-07-20
### Added
- genome input for bonus
- option to to the 'percent' method for determining taxonomy
- testing for taxonomy
- Tdb now includes the columns "best_hit" and "full_tax" (for both methods of determining taxonomy)

## [1.0.0] - 2017-07-14
### Added
- bonus --check_dependencies now exists
- gANI now gives a message when it fails

### Fixed
- pytest is now automatically installed with pip
- the logger now does ' '.join(args) when printing the args to run dRep
- documentation now correctly says "ANIcalculator"
- makes sure ANIcalculator and checkM work when loading them (even if you can find them in the system path)

## [0.5.7] - 2017-04-26
### Fixed
- having the user add their own Chdb now works
- in bonus, changed an erroneous > to >= when calling centrifuge
- changed the way the test_suite compares dataframes

## [0.5.6] - 2017-04-23
### Fixed
- nc option now works with gANI (controller makes it a float)
- test_backend is now a thing automatically

## [0.5.5] - 2017-04-09
### Fixed
- coverage values at the threshold are now accepted (< instead of <=)
- setup now automatically installs scikit-bio (needed for MDS plot)

### Changed
- the loop in which genome lengths are calculated in d_cluster is changed to
prevent errors when running large numbers of genomes

## [0.5.4] - 2017-03-28
### Changed
- default coverage method is now larger (tests are updated to reflect this)

### Fixed
- gANI now properly computes the coverage using the "larger" method
- gANI can now tell when it's not installed

## [0.5.3] - 2017-03-28
### Changed
- documentation about dependencies changed; versions added as well
- changed the way ani averaging is done; substantial speed increase when working with very large secondary clusters
- added the dRep version to the log file

## [0.5.2] - 2017-03-24
### Fixed
- prodigal now properly threads for gANI

## [0.5.1] - 2017-03-17
### Fixed
- Plot 3 (MDS) now prints onto a square grid
- Testing suite now launches all tests automatically

## [0.5.0] - 2017-03-16
### Added
- The final analyze plot can now be made!

### Fixed
- pyenv message
- log typos

## [0.4.0] - 2017-03-10
### Fixed
- prodigal now multithreads the correct amount
- compare_wf double-printing issue
- --SkipMash now works

### Changed
- default checkM method is now lineage_wf in dereplicate_wf, filter, and choose
- pyenv now reverences non-anaconda in the documentation and error messages
- default min length in now 50,000

### Added
- more thorough testing (though still not nearly enough)
- log now contains the exact command run

## [0.3.8] - 2017-03-06
### Fixed
- Mash sketch size argument actually works now

## [0.3.7]
### Fixed
- Fixed analyze to produce non-ugly plots

## [0.3.6]
### Added
- Can now change the method of calculating coverage with ANIn
- Added TODO

## [0.3.5]
### Fixed
- Bug that caused the program to crash when looking for Mash
- Lots of typos (in documentation and program help)

## [0.3.4]
### Changed
- Fixed some typos in documentation
- Added bioRxiv link in documentation

### Fixed
- Edited matplotlib to make text that can be edited by illustrator

## [0.3.3]
### Fixed
- Small bug in telling the user when checkM/prodigal is not in system path
- Small bug preventing SkipSecondary from working

### Added
- Test for SkipSecondary

## [0.3.2]
### Fixed
- Time estimates now take into account threading
- Fixed the "choose" operation erroneously taking the "adjust" arguments as well
- Removed the 'mauve' option as a secondary algorithm
- Fixed a lot of erroneous logging

### Changed
- Removed a bunch of commented-out methods

## [0.3.1]
### Changed
- READme updated

### Added
- Documentation

### Fixed
- A bug where '-h' wasn't working properly
- Some messages in cluster that should be going to log

## [0.3.0] - 2017-01-30
### Added
- Basic functional testing
- This here Changelog
- Probably a lot of other stuff that I don't remember (due to the previous lack of Changelog)

### Changed
- Moved the argument parsing to a "controller" class
- The backend of the way the log runs. It's better now

## [0.1.0] - 2017-01-15
- This is the name of the release generated before implementation of a changelog. The development version number for this release was 0.2.3
