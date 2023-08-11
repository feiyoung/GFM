# GFM

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/GFM)](https://cran.r-project.org/package=GFM)
[![](https://cranlogs.r-pkg.org/badges/GFM?color=orange)](https://cran.r-project.org/package=GFM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/GFM?color=orange)](https://cran.r-project.org/package=GFM)
<!-- badges: end -->

GFM: Generalized factor model for ultra-high dimensional variables with mixed types.

GFM  is a package for analyzing  the (ultra)high dimensional data with mixed-type variables, developed by the Huazhen Lin's lab. It is not only computationally efficient and scalable to the sample size increment, but also is capable of choosing the number of factors. In our JASA paper, a two-step method is proposed to estimate the factor and loading matrix, in which  the first step used the alternate maximization (AM) algorithm to obtain initial estimator. In the paper, the information criterion was provided to determine the number of factors.  Recently, we proposed an overdispersed generalized factor model (OverGFM) and designed a variational EM algorithm to implement OverGFM. A  singular value ratio based method was provided to determine the number of factors. In addition, the estimate from OverGFM can be also used as the initial estimates in the first step for GFMs in our previous JASA paper. 


Check out our [JASA paper](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999818?journalCode=uasa20) for alternate maximization and information criterion, and our [Package vignette](https://feiyoung.github.io/GFM/docs/index.html)  for a more complete description of the usage of  GFM  and OverGFM. 

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


# NEWs
GFM version 1.2.1 (2023-08-10)

The function `overdispersedGFM()` that implements the overdispersed generalized factor model is added. 
In addition, the function `OverGFMchooseFacNumber()` is added, which implements singular value ratio (SVR) based method to select the number of factors.


