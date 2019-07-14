# gemma2 0.0.2

## Major changes

* None

## Minor changes

* New version number for CRAN submission



# gemma2 0.0.1.1

* added an input parameter, verbose_output, to MphEM function. The default value is verbose_output = FALSE, which affects the output by returning only the last iteration's parameter values. In fact, it does more than that by actually overwriting each iteration's parameter values with the subsequent iteration's values. To recover the 'old' behavior, use verbose_output = TRUE. I also added new tests to ensure that the value of the verbose_output argument doesn't affect the last iteration's values. These tests are in tests/testthat/test_MphEM.R.

# gemma2 0.0.1

* added vignette for comparing mvlmm gemma2 results with those of Zhou's C++ GEMMA code (v 0.96)

# gemma2 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.



