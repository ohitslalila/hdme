# hdme 0.1.1.9002
Added a `plot.gds` function, which plots the coefficients estimated by `fit_gds`.

# hdme 0.1.1.9001
Internal adjustment, which removed importing of external packages into the namespace, and instead specified the functions explicitly using `::`.

# hdme 0.1.1
`tidyverse` has been removed from **Suggests** field `DESCRIPTION`. `dplyr` and `tidyr` have been added instead. Similarly, `library(tidyverse)` in the vignette has been replaced by `library(dplyr)`, `library(tidyr)`, and `library(ggplot2)`.

# hdme 0.1.0
`hdme` is now on CRAN.