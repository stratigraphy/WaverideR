# Changes in version 0.3.0

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



# WaverideR 0.2.0

* Added a `NEWS.md` file to track changes to the package.
