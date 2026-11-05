# An ETP astronomical solution

The `astrosignal_example` is a pre-generated ETP  
(eccentricity-tilt-precession) (p-0.5t based on the la2004 solution)  
the `astrosignal_example` can be used to anchor the
[`grey`](https://stratigraphy.github.io/WaverideR/reference/grey.md)  
data set to an astronomical solution eg. `astrosignal_example`  
using the
[`astro_anchor`](https://stratigraphy.github.io/WaverideR/reference/astro_anchor.md)
function. the data set was generated using the  
[etp](https://rdrr.io/pkg/astrochron/man/etp.html) function of the
'astrochron' R package. The pre-generated ETP spans 5000 to 6000kyr.

## Details

Column 1: time (kyr)  
Column 2: ETP  

## References

Stephen R. Meyers,Cyclostratigraphy and the problem of astrochronologic
testing, Earth-Science Reviews,Volume 190,2019,Pages 190-223,ISSN
0012-8252 \<doi:10.1016/j.earscirev.2018.11.015\>

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B.
Levrard, B., 2004, A long term numerical solution for the insolation
quantities of the Earth: Astron. Astrophys., Volume 428, 261-285.
\<doi:10.1051/0004-6361:20041335\>

## Author

Generated using the [etp](https://rdrr.io/pkg/astrochron/man/etp.html)
function of the
[astrochron-package](https://rdrr.io/pkg/astrochron/man/astrochron-package.html).
