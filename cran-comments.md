## kdml -- Resubmission

This is a resubmission. The following four notes were sent by CRAN on August 21st 2024.

# 1st comment

Please write references in the description of the DESCRIPTION file in
the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: authors (year) <[https:...]https:...>
with no space after 'doi:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")
So please add some form of linking as well.

We have fixed the two erroneous references to adhere to the first form.  

# 2nd comment

Please omit the <<tt>> <<\tt>> tags as they are not rendered properly.
Please add () behind all function names in the description texts
(DESCRIPTION file). e.g: --> dkss(), dkps()

We have fixed the two function names in the description to the correct form.

# 3rd comment

You write information messages to the console that cannot be easily
suppressed.
It is more R like to generate objects that can be used to extract the
information a user is interested in, and then print() that object.
Instead of print()/cat() rather use message()/warning() or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console. (except for print, summary, interactive functions)
-> R/dkps.R; R/dkss.R; R/mscv.dkps.R; R/mscv.dkss.R

# 4th comment

Almost all your examples are wrapped in \donttest{} and therefore do not
get tested.
Please unwrap the examples if that is feasible and if they can be
executed in < 5 sec for each Rd file or create additionally small toy
examples to allow automatic testing.

We removed the \donttest{} and compiled the package on the test environments and 
received the following note:

❯ checking examples ... [110s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
  mscv.dkps 24.63   3.38   28.08
  dkss      24.80   2.31   27.36
  dkps      24.92   2.15   27.64
  mscv.dkss 21.05   2.67   23.77
  
These times are > 5 sec, and thus we kept the \donttest{} in the manual file.

## Test environments

- Local: Windows 10, R 4.4.1 (x86_64-w64-mingw32/x64)
- win-builder: Windows Server 2022 x64 (build 20348, x86_64-w64-mingw32)
- GHA Builder:
  - macOS Sonoma 14.6-R version 4.4.1 (2024-06-14)
  - Ubuntu 22.04.4 LTS-R version 4.5.0
  - Ubuntu 22.04.4 LTS-R version 4.3.3 (2024-02-29)
  - Ubuntu 22.04.4 LTS-R version 4.4.1 (2024-06-14)
  - Windows Server 2022 x64 (build 20348)-R version 4.4.1 (2024-06-14)

## R CMD check results (Aug 21)

0 errors | 0 warnings | 2 notes

❯ checking CRAN incoming feasibility ... [23s] NOTE
  Maintainer: 'John R. J. Thompson <john.thompson@ubc.ca>'
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Ghashti (7:472, 7:562)

This is not a misspelling.

❯ checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

This is not a problem associated with this package. 