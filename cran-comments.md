# Changes in version 0.3.3

fixed and error in: track_period_wavelet


# Changes in version 0.3.2

since submission:
## CRAN package check 
Kurt Hornik writes:

Still open ...

Please fix before 2023-09-10 to safely retain your package on CRAN.

Dear maintainer,

R-devel will soon find

Found the following (possibly) invalid DOIs:
  DOI: 1016/j.epsl.2018.02.010
    From: man/mag.Rd
    Status: 404
    Message: Not Found
  DOI: <doi:10.1190/1.1441328>
    From: man/analyze_wavelet.Rd
          man/plot_wavelet.Rd
    Status: 404
    Message: Not Found
  DOI: <doi:10.1190/1.1441329>
    From: man/analyze_wavelet.Rd
    Status: 404
    Message: Not Found
  DOI: doi.org/10.1051/0004-6361/201117504
    From: man/astro_anchor.Rd
    Status: 404
    Message: Not Found

for you package: previously DOIs in Rd files were missed in the CRAN
submission checks.

Please fix.

When doing so, please note that for the valid DOIs you do things like

#win_fft.Rd:time series analysis <\doi{doi:10.1016/j.earscirev.2018.11.015}>

#where you enclose the \doi{} inside <...>: please don't do this (and
#drop the extra doi:).  One needs to do e.g.

 <doi:10.1016/j.earscirev.2018.11.015>

 in the package DESCRIPTION Description, but should do

 \doi{10.1016/j.earscirev.2018.11.015}

in Rd.

Please fix before 2023-08-23.

Best
-k

The (possibly) invalid links were fixed
The way links should be written in the documentation and documentation
is also adjusted to be inline with CRAN guidelines


added functions:

add_wavelet
adds a wavelet plot to a larger plot

add_wavelet_avg
adds the average wavelet power to a larger plot

curve2time_unc
converts a tracked period curve to time taking the uncertainty of the curve 
into account

curve2time_unc_anchor
converts a tracked period curve to time by taking the uncertainty of the curve 
and anchored radiometric date into account 

dur_gaps
calculates the duration off a hiatus based on the missing part of an astronomical
cycle

geo_col
generate standard color codes for the Geological Time Scale

geo_loc
generates ages for the boundaries of a geochronological subdivision

geo_mid
generate the mean age of a geological subdivision

lithlog_disc
discriticizes lithologs from a block pattern to a continous signal which
allows it to analyzed using spectral analysis techniques

minimal_tuning
Creates an age model using the minimal tuning technique.
This means that the distance between 2 peaks of an extracted cycle is 
set to duration of the interpreted astronomical cycle

retrack_wt_MC
Re-track cycles using a Monte-Carlo simulation
Then input is an already tracked period of an astronomical cycle and the wavelet
in which it was tracked. the already tracked periods are combined to create a 
pseudo-curve which is then re-tracked in one of the wavelets (randomly selected)
the result is a bunch of retracked curves 

sedrate2tune
converts a data set using the sedimentation rate to a tuned age-model

updated functions: 

Hilbert_transform
added side padding to avoid edge effects

plot_wavelet
added additional color plotting option 

plot_win_fft
added additional color plotting option 

track_period_wavelet
added additional color plotting option 
added the option to select any point within the plot

#Changes in version 0.3.1

Examples contained dontrun{} which were replaced with donttest{} when examples are executed in >5 seconds
otherwise dontrun{} is removed completely. Searched all files using CTR:+SHIFT+F to see if dontrun{} was still
present and dontrun{} was not found after replacement

Checked and can now confirm that only 1 core is used by default in functions which can run
using multiple cores by default the parameter to run multiple cores is set to (run_multicore = FALSE).
In examples run_multicore is set to; (run_multicore = FALSE).

Information messages were still written to the console which is replaced with the if(verbose)cat(..) option
for interactive scripts where printing text to the console is important text is added to the examples and arguments of text indicating that functionality is reduced if verbose is set to FALSE meaning that only basic functionality is
maintained

package names, software names and API are now written in quotes e.g. '----'.

some plotting functions opened new windows for new plots, these plot would not show up in Rmarkdown documents.
An option was added to turn off the open plots in new windows. 

added support for more color palette options.
added support to switch the x and y-axis when plotting.

## CRAN submission feedback to version 0.3.0

Please always write package names, software names and API (application 
programming interface) names in single quotes in title and description. 
e.g: --> 'WaverideR'
Please note that package names are case sensitive.

You still have example wrapped within dontrun{} which I think is not 
necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace 
dontrun{} with donttest{} or explain why dontrun{} is necessary for 
these examples.

You still write information messages to the console that cannot be 
easily suppressed. It is more R like to generate objects that can be 
used to extract the information a user is interested in, and then 
print() that object.
Instead of print()/cat() rather use message()/warning()  or 
if(verbose)cat(..) (or maybe stop()) if you really have to write text to 
the console.
(except for print, summary, interactive functions)

Please confirm that  you do not use more than 2 cores in your examples, 
vignettes, etc.

##Changes in version 0.3.0

General changes: 

* Added References both in the description file of the package and the description texts of individual functions.

* Changed print() statements to if(verbose)cat(..) to allow suppression of messages in the console
An exception was made for the astro_anchor function where cat() is used instead because because "feedback" needs
to entered into the console

* spell-checked code and removed misspelled words 

* Removed dontrun from examples and replaced with donttest when examples are executed in >5 seconds
otherwise dontrun is removed and extra code is added. Many examples still have donttest added to them as the examples of functions use real-world examples of already published records. Using real-world example data-sets is preferred over synthetic data sets as  

* Previous iterations of functions often contained code which changed the user's options. by adding oldpar <- par(no.readonly = TRUE) followed by  on.exit(par(oldpar)) to functions this issue is avoided. 
This however makes images/plots uneditable after plotting. To keep images/plot editable the parameter "keep_editable" is added to functions which generate images/plots, If keep_editable is set to TRUE plot/images remain editable but by default keep_editable is set to FALSE to avoid changing the user's options.

* Functions which  have the option to run on  multiple cores have the parameter "run_multicore" added to them. On default run_multicore is set to FALSE to avoid running the process on multiple cores overloading the user's computer. In all examples run_multicore is set to FALSE meaning that only 1 core is used when a function is executed. 

Function specific changes:

* added astrosignal_example data set
* added depth_rank_example data set 
* added extract_amplitude function 
* added flmw function function
* added lithlog_disc function
* added plot_sed_model function
* added plot_win_fft function
* added win_fft  function
* major revision of the max_detect and min_detect functions which reduced the possibility of generating erroneous output

## CRAN submission feedback

If there are references describing the methods in your package, please 
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for 
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

Possibly misspelled words in DESCRIPTION:
         cyclcostratigraphic (11:137)
-> probably it should be cyclostratigraphic

\dontrun{} should only be used if the example really cannot be executed 
(e.g. because of missing additional software, missing API keys, ...) by 
the user. That's why wrapping examples in \dontrun{} adds the comment 
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace 
\dontrun{} with \donttest{}.

You write information messages to the console that cannot be easily 
suppressed. It is more R like to generate objects that can be used to 
extract the information a user is interested in, and then print() that 
object.
Instead of print()/cat() rather use message()/warning()  or 
if(verbose)cat(..) (or maybe stop()) if you really have to write text to 
the console.
(except for print, summary, interactive functions)

Please make sure that you do not change the user's options, par or 
working directory. If you really have to do so within functions, please 
ensure with an *immediate* call of on.exit() that the settings are reset 
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
If you're not familiar with the function, please check ?on.exit. This 
function makes it possible to restore options before exiting a function 
even if the function breaks. Therefore it needs to be called immediately 
after the option change within a function.

Please ensure that you do not use more than 2 cores in your examples, 
vignettes, etc.

Please fix and resubmit.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

#

