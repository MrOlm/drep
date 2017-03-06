# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [UNCOMMITTED]

## [0.3.8]
### Fixed
- Mash sketch size argument actually works now

## [0.3.7]
### Fixed
- Fixed anaylze to produce non-ugly plots

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
