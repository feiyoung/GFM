# GFM

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/GFM)](https://cran.r-project.org/package=GFM)
[![](https://cranlogs.r-pkg.org/badges/GFM?color=orange)](https://cran.r-project.org/package=GFM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/GFM?color=orange)](https://cran.r-project.org/package=GFM)
<!-- badges: end -->

GFM: Generalized factor model for ultra-high dimensional variables with mixed types.

GFM  is a package for analyzing  the (ultra)high dimensional data with mixed-type variables, developed by the Huazhen Lin's lab. It is not only computationally efficient and scalable to the sample size increment, but also is capable of choosing the number of factors. Two algorithms, variational EM and alternate maximization, are designed to implement the generalized factor model, respectively, which ensures that the factor matrix and loading matrix together with the number of factors can be well estimated. 


Check out our [JASA paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999818?journalCode=uasa20) and our [Package vignette](https://feiyoung.github.io/GFM/docs/index.html)  for a more complete description of the methods and analyses. 

GFM can be used to analyze experimental dataset from different areas, for instance:

* Social and behavioral sciences
* Economy and finance
* Genomics...


Please see our new paper for model details:

[Wei Liu, Huazhen Lin, Shurong Zheng & Jin Liu (2021) . Generalized factor model for ultra-high dimensional mixed data. Journal of the American Statistics Association (Online).](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999818?journalCode=uasa20)

# Installation

To install the the packages 'GFM' from 'Github', firstly, install the 'remotes' package.
```{Rmd}
install.packages("remotes")
remotes::install_github("feiyoung/GFM")
```
Or install the the packages "GFM" from 'CRAN'
```{Rmd}
install.packages("GFM")
```

# Demonstration

For an example of typical GFM usage, please see our [Package vignette](https://feiyoung.github.io/GFM/docs/index.html) for a demonstration and overview of the functions included in GFM.
