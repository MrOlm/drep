# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [4.0.0] - Unreleased

dRep v4 makes genome clustering scale. The headline change is that primary and
secondary clustering now both default to skani, and run as a **single pass** over
the data instead of comparing every genome twice.

**This release changes results.** Read "Breaking changes" below before upgrading
an existing analysis. If you need the old behavior, `--primary_algorithm MASH
--S_algorithm fastANI --primary_clusterAlg average` gets close, but the skani
coverage fix (see Fixed) cannot be turned off, and it was a genuine bug.

### Breaking changes

- **skani is now the default for both clustering steps** (`--primary_algorithm`
  defaults to `skani`, was MASH; `--S_algorithm` defaults to `skani`, was
  fastANI). skani must be installed. Mash is now only needed for
  `--primary_algorithm MASH`, and the startup dependency check no longer demands
  it otherwise.
- **Primary clustering now defaults to single linkage.** Previously `--clusterAlg`
  (default `average`) drove *both* clustering steps. Primary now has its own
  `--primary_clusterAlg`, defaulting to `single`. Single linkage is the right
  choice for a deliberately inclusive pre-filter, and it is what makes the
  low-memory algorithm possible. `--clusterAlg` still controls secondary
  clustering and still defaults to `average`.
- **`Mdb.csv` means something different under `--primary_algorithm skani`.** It is
  now a *sparse* table of real skani ANI values plus alignment coverage, holding
  only pairs above skani's screening threshold — roughly 779k rows for 10,000
  genomes, versus 100M rows of dense Mash distances. Genomes with no
  above-threshold pairs do not appear in it at all. Anything parsing `Mdb.csv`
  needs to account for this.
- `--S_algorithm skani` results change; see the coverage fix under Fixed.
- `--low_ram_primary_clustering` was removed (see Removed).

### Added

- `--primary_algorithm {skani,MASH}` — choose the primary clustering program.
- `--primary_clusterAlg` — linkage method for primary clustering, independent of
  the secondary `--clusterAlg`.
- `--classic_primary_clustering` — force the pre-v4 dense scipy primary path.
- `--primary_skani_min_af` — minimum aligned fraction for a pair to form a
  primary-clustering edge (skani only, default 15).
- `--no_reuse_primary_comparisons` — re-run skani during secondary clustering
  rather than reusing primary's comparisons. A debugging escape hatch; reuse is
  exact.

### Changed

- **Primary clustering no longer builds the N x N distance matrix.** Single-linkage
  clustering at a fixed threshold is identical to finding connected components,
  so it is now computed directly with union-find, removing the pivot-then-unpivot
  RAM spike (issue #259). Memory is O(genomes + edges) instead of O(genomes^2).
  On a synthetic 8,000-genome set, peak RAM for this step dropped from 7.42 GB to
  0.55 GB; the old path grew quadratically while the new one stays flat. Together
  with `--primary_algorithm skani`, this addresses the out-of-memory crashes
  reported when clustering tens of thousands of genomes.
- **Secondary clustering reuses primary's comparisons.** With skani for both
  steps, dRep previously sketched every genome twice and computed the same ANI
  values twice: once across all genomes, then again within each primary cluster.
  Since secondary only compares genomes *within* a primary cluster, those pairs
  are a subset of what primary already computed. On 10,000 UHGG genomes, 94% of
  the pairs driving secondary clustering were already present with identical ANI
  to 6 decimal places. dRep now runs skani once and derives both steps from it.
  Secondary clustering went from 15 minutes to 24 seconds, producing an identical
  partition; whole-pipeline `dereplicate` went from 22.4 to 13.8 minutes.
- Primary clustering with skani writes far less to disk: 82 MB vs 15 GB of Mash
  output for 10,000 genomes.
- `--multiround_primary_clustering` and `--primary_chunksize` now warn that they
  only apply to `--primary_algorithm MASH`. skani's sparse output never builds
  the N x N table that multiround exists to avoid, and has none of multiround's
  chunk-splitting imprecision.
- `--SkipMash` help text clarified: it skips primary clustering whatever the
  primary algorithm is. The name is historical.

### Removed

- `--low_ram_primary_clustering`. Union-find is now the default for single-linkage
  primary clustering, so the flag had become a no-op.
- **networkx is no longer a dependency.** It was only used by the connected-components
  path behind `--low_ram_primary_clustering`.

### Fixed

- **skani alignment coverage was a percent, not a fraction (results-affecting).**
  skani reports aligned fractions as 0-100, but `load_skani` only divided ANI by
  100 and passed the aligned fraction through untouched. Every other algorithm
  reports `alignment_coverage` on a 0-1 scale, which is the scale `cov_thresh` is
  compared against, so the coverage filter was effectively inert for
  `--S_algorithm skani`: a pair aligning over 1% of the genome had
  `alignment_coverage=1.04` and sailed past a `cov_thresh` of 0.5. On the bundled
  test genomes this merged *E. casseliflavus* with *E. faecalis* — two different
  species — into one secondary cluster. Coverage filtering now actually applies.
- `ScaffoldLevel_dRep.py` crashed on MUMmer 3 with "nucmer failed with exit code
  1". It passed `-t` (threads) unconditionally, but that option only exists in
  MUMmer 4, and MUMmer 3 rejects it rather than ignoring it. `conda install
  mummer` still installs 3.23. The script now detects whether nucmer supports
  `-t` and only passes it if so; MUMmer 4 keeps its threading.
- The primary dendrogram is still produced for modest genome sets under the
  streaming/sparse paths, which build no linkage matrix of their own. Above
  `--primary_dendrogram_max_genomes` (2000) it is skipped, as multiround already
  did.

### Validation

v4 clustering was validated against 10,000 real genomes from the UHGG catalogue
(24 GB, 1,248 species, 9,599 MAGs + 401 isolates). dRep independently recovered
1,232 secondary clusters at `-sa 0.95`, against UHGG's own 1,248 species
assignments: 2.4% of clusters spanned more than one UHGG species, and 2.9% of
UHGG species were split across clusters. The one-pass path reproduced the
two-pass result exactly — same primary clusters, same secondary clusters, same
representative genomes.

## [3.7.1] - 2026-06-30
- Fix crash when fewer than 2 genomes remain after filtering (issue #300)
- Fix argument list bug (issue #288)
- Dendrogram label contrast fix (issue #235)
- Fix Widb.csv showing NA for completeness/contamination when using --genomeInfo (issue #220)
- Clarify --genomeInfo help text: genome name must include file extension (issue #220)
- Fix plotting crash with tertiary clustering due to decimal cluster names (issue #202)
- Add warning when MASH sketch size may be too small for low P_ani thresholds (issue #137)
- Fix false "fastANI isn't working" error when fastANI -h returns non-zero exit code (issue #96)

## [3.7.0] - 2026-06-30
- Implement testing with the GitHub CLI
- Maintain compatibility with pandas 3

## [3.6.2] - 2024-04-14
- update requirements

## [3.6.1] - 2024-03-19
- Add low_ram_primary_clustering option for memory-efficient primary clustering
- Add tests for low_ram_primary_clustering functionality
- Update documentation to include new low_ram_primary_clustering option

## [3.5.1] - 2024-03-11
- Update goANI to only include 1 best query 

## [3.5.0] - 2024-02-21
- Lots of little bug-fixes for scoring plots to update to new versions of python / matplotlib
- Remove Biopython as a dependency
- Add tqdm as a dependency in pip
- Add skani as a comparison option

## [3.4.5] - 2023-09-15
- Tiny bugfix related to loading cached secondary clusters

## [3.4.4] - 2023-09-14
- Tiny bugfix related to loading cached secondary clusters

## [3.4.3] - 2023-04-18
- Fix pandas pivot bug 

## [3.4.2] - 2023-02-06
- Quick bugfix for last update 

## [3.4.1] - 2023-01-11
- Only run checkM on bugs not in the genomeInfo list
- https://github.com/MrOlm/drep/pull/179
- Thanks Jon Sanders!

## [3.4.0] - 2022-07-31
- Add the options "--skip_evaluate" and "--skip_analyze" to dereplicate

## [3.3.1] - 2022-07-15
- Update parse_stb.py to allow a list of genomes as the .fasta input in --reverse mode

## [3.3.0] - 2022-06-13
- Set default S_algorithm to fastANI and sa to 0.95
- Run Prodigal in "single" mode
- If using multiround_primary_clustering, yell if not also tertiary_clustering

## [3.2.2] - 2021-06-22
- Better timing prediction when using fastANI

## [3.2.1] - 2021-06-14
- clusterAlg is now used to cluster the initial groups when doing mutliround_primary_clustering (before it was always "single")
- fixed a bug that happened when you ran skipSecondary with centW > 0 (https://github.com/MrOlm/drep/issues/120)

## [3.2.0] - 2021-03-09
- checkM can now be run in groups (added argument variable "checkm_group_size")
- other slight internal restructuring around how checkM works
- a little more info added to the "troubleshooting dRep" section

## [3.1.1] - 2021-03-04
- Update parse_stb.py to "append" instead of open, which lets it handle way bigger numbers of genomes
- Report how many genomes are being compared in last step of multi-round primary clustering

## [3.1.0] - 2021-02-25
- Add the "extra_weight_table" flag to allow uses to add custom extra scores to their genomes
- Do some automatic checking of genome input (List number of genomes, throw warning if same base name is used twice, throw warning if over 5000 genomes and no —multiround_primary_clustering)
- Make a section in the docs for troubleshooting checkM and reference it when checkM fails

## [3.0.1] - 2021-01-27
- Remove some assert statements
- Really just a version bump to update bioconda requirements

## [3.0.0] - 2021-01-07
- Refactoring the test suite and the d_cluster module
- Adding help to the -g option
- Bare-bones support for gzipped genomes (lots of dependencies dont handle them)
- Give a warning when it has errors loading the Mash table
- Completely refactor the test suite to use pytest
- Make plotting only give tracebacks when run in debug mode
- Remove most of the options (just `dereplicate` and `compare` remain)
- Add greedy clustering support! Both `multiround_primary_clustering` and `greedy_secondary_clustering`
- Add centrality support; also handle centrality with greedy_secondary_clustering (will be calculated with Mash)
- Add --run_tertiary_clustering. This feature runs a final dRep job within the original dRep folder,
  and adjusts Cdb and Wdb accordingly (see `run_tertiary_clustering` in `d_evaluate.py`)

## [2.6.2] - 2020-05-28
- Log information about GenomeInformation when loading it

## [2.6.1] - 2020-05-22
- Numerous improvements to ScaffoldLevel_dRep.py, including ability to process in chunks
- Update parse_stb.py to handle zipped .fasta files

## [2.6.0] - 2020-05-14
- Add helper scripts ScaffoldLevel_dRep.py and parse_stb.py

## [2.5.4] - 2020-03-06
- Trying to fix a bug related to pandas categories

## [2.5.3] - 2020-02-27
### Fixed
- More bug fixes related to FastANI
- Allow loading of cached Ndb.csv

## [2.5.2] - 2020-02-26
### Fixed
- The bug I tried to fix in 2.5.1 is able to fix itself

## [2.5.1] - 2020-02-24
### Fixed
- Instead of crashing out, FastANI will report the error and keep going if parsing fails

## [2.5.0] - 2020-02-20
### Added
- FastANI is now an option for secondary comparisons
- You can now feed in a list of genomes via the -g option

## [2.4.2] - 2020-02-01
### Changed
- More edits to make goANI work (what a bizarre bug with the output sometimes having different headers?)

## [2.4.1] - 2020-01-23
### Changed
- Changed the flag -n_PRESET to --n_PRESET
- Handle the case where a nsimscan run completely fails in goANI mode
- Remove "--force_overwrite" from checkM since its no longer supported

## [2.4.0] - 2020-01-07
### Changed
- Updated warnings and documentation to reflect checkM being in python3 now
- Updated help to link to documentation
- Updated documenation in other little ways

## [2.3.2] - 2019-03-28
### Fixed
- print tracebacks when plots fail
- fixed a weird bug with plot 5 resulting from genomeInformation having too many columns
- fixed plot 6 failing due to the deprecatoin of d_cluster.av_ani

## [2.3.1] - 2019-03-09
### Fixed
- goANI bug resulting when there is no overlap in a filtered nsimscan

## [2.3.0] - 2019-02-18
### Added
- goANI is now added as a secondary clustering algorithm; an open-source alternative to gANI

## [2.2.4] - 2019-01-31
### Changed
- Renamed --noQualityFiltering to --ignoreGenomeQuality
- Changed some things around to satisfy pandas deprication warnings

## [2.2.3] - 2018-10-11
### Fixed
- Fixed bug where Mash dendrogram labels were scrambled if a big list of genomes was used (thanks brymerr921)
- Fixed typos (thanks AstrobioMike!)

## [2.2.2] - 2018-06-18
### Added
- Added the --set_recursion option for filter to handle dendropy errors

## [2.2.1] - 2018-05-30
### Changed
- WorkDirectory now loads databases in a way that makes more sense for large tables
- Some extra caching debug options
- Some commented out memory stuff
- Mash comparisons are now actually multithreaded (thanks mruehlemann)
- Throws an error if run with python2

## [2.2.0] - 2018-04-02
### Changed
- RAM optimization with regards to loading Mash table
- Pickle protocol 4 to allow storage of larger clusterings
- Increased debug tools

## [2.1.1] - 2018-03-22
### Changed
- use threading instead of multiprocessing. This should significantly help with RAM utilization of large genome lists
- Added some extra debug options

## [2.1.0] - 2018-03-19
### Changed
- removed the overwrite option and enabled it by default. It was half-baked and didn't work anyways

## [2.0.6] - 2018-03-16
### Added
- added unit tests for scoring

## [2.0.5] - 2018-01-12
### Fixed
- taxonomy now works when prodigal was not previously run

## [2.0.4] - 2018-01-09
### Fixed
- some changes in d_cluster that make gANI work
- changes to the test suite so it doesn't fail if centrifuge isn't installed

## [2.0.3] - 2018-01-08
### Fixed
- Fixed test_suite to work

## [2.0.2] - 2018-01-04
### Fixed
- Fixed try / excepts around the plots failing. Looks like a matplotlib issue

## [2.0.1] - 2018-01-02
### Fixed
- Plot 6 will not crash if Widb is not present
- Put try / excepts around the plots failing

## [2.0.0] - 2017-12-10
### Added (general)
- API documentation

### Changed (general)
- dereplicate_wf to dereplicate
- compare_wf to compare
- removed adjust (simply removed reference to it in the argument parser; easy to bring back if desired)

### Added to d_filter
- Complete API coverage
- Ability to include genome information from other sources more easily
- More tests

### Added to d_cluster
- argument parsing broken up into more groups
- more tests
- complete API coverage
- ANImf, gANI, and ANIn now make folders in the output, so that it doesnt have too
many files in one dir
- mash paste now goes in chunks, so that it will work if you have huge numbers of genomes (getting around OSError: [Errno 7] Argument list too long)

### Added to d_choose
- full API coverage
- ability to make / generate genomeInfo better

### Added to d_analyze
- full API coverage
- deleted the whole re-cluster thing

## [1.4.3] - 2017-10-24
### fixed
- fixed centrifuge call (shell=False)

## [1.4.2] - 2017-10-24
### Changed
- made ANImf the default comparison algorithm
- added a little documentation on what ANImf is

### fixed
- made ANImf not rerun if .delta.filtered is present
- made bonus call centrifuge using the new calling method
- made the new calling method actually take into account the number of threads
- updated tests to account for default ANImf

## [1.4.1] - 2017-10-23
### Fixed
- with gANI, fixed a minor bug with the naming of self-comparisons in Ndb.csv
- with ANImf, fixed a major bug that was preventing it from working right at all

## [1.4.0] - 2017-10-20
### Added
- added ANImf comparison method

## [1.3.2] - 2017-10-20
### Fixed
- test_suite now works
- removed unnecessary import

## [1.3.1] - 2017-10-19
### Changed
- output of all external commands is now stored in the log directory
- all external commands are now run through the same method

### Fixed
- mash should now work on larger argument lists

## [1.3.0] - 2017-10-05
### Fixed
- filtering with regards to strain heterogeneity fixed

## [1.2.1] - 2017-10-03
### Changed
- Logging now logs mummer commands run with ANIn option

### Added
- Choose can now function without checkM

## [1.2.0] - 2017-09-11
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
