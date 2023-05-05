#Changes in version 0.3.0

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

# CRAN submission feedback

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

