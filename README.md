# Orientations

R code for doing fabric analysis.  Initially intended for use with total station coordinate data (two XYZ points on the long axis of an object) but bearing and plunge data can be transformed into the appropriate format for use here (see below).

This repository contains the code and markdown document (with supporting files) previously published in PLOS One (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190195).  Some small bug fixes have been added (see below).  The rMarkdown document (McPherron 2018 - PLOS One.Rmd) is the paper itself and when knitted in R will produce a PDF of the paper.  This document shows how the set of functions found in orientations.R can be used to do fabric analysis.  Otherwise, currently there is no manual for the functions found in orientations.R.

If you find problems with this code, I welcome bug fixes.  Please either contact me or if you are familiar with GitHub you can use a pull request.  If the changes are substantial, I would appreciate being contacted first.

Please cite *McPherron 2018* if the code is used in publications.

### Converting bearing and plunge angles for use with this code

The orientation.R code assumes that the data come from a total station with two XYZ coordinates for each object.  However, as many researchers have data as bearing and plunge angles only, I have written the code in the Convert Bearing and Plunge folder to convert these data into a format that orientations.R can handle.  For now I cannot provide the data that I used to test this code, but when it becomes available I will.

Please, if you find an errors with this approach, let me know.

### Bug/Feature Fixes

- In the original PLOS One paper I had added a small feature to highlight objects with high isotropy.  I have since removed this (thanks to Sam Lin for spotting this feature, I mean bug).
- I changed package references from require() to library() to make it easier to know what libraries to install.