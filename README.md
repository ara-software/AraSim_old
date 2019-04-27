# AraSim_old

## Warning
This is a *legacy* version of AraSim. This repository will be archived and will not support new issues, branches, pull requests, etc. Instead, you should do work on the active development repository of AraSim: <https://github.com/ara-software/AraSim>.

## Why an "old" version?
Back when AraSim was hosted on an SVN at Ohio State, in March of 2014, it underwent a directory structure change. Prior to March 2014, all commits were made directly to a monolithic "central" folder in the repository.

In March 2014, Carl Pfender changed from this monolithic structure to the `trunk`, `branches`, and `releases` structure that is more typical of SVN repositories. Because of this, it is difficult to connect the commit history of the repository when it was monolithic to the commit history when it was under the more standard structure.

So, we settled on having two repositories. This `AraSim_old` version captures r1-446, and the normal `AraSim` version is r446-803 (803 was the last SVN revision before the migration to GitHub). `AraSim_old` is used purely for historical records; all modern development will happen in `AraSim`.