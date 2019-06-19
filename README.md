# slm
---
Overview of the slm package
---

The slm package enables to fit linear models on datasets considering the dependence between the observations.
Most of the functions are based on the functions and methods of lm, with the same arguments and the same format for the outputs.

---
The main functions
---

* slm function, in "slm-main.R":
The slm function is the main function of this package. Its architecture is the same as the lm function
but it takes into account the possible correlation between the observations. To estimate the asymptotic covariance matrix of
the least squares estimator, several approaches are available: "fitAR" calls the
cov_AR function, "spectralproj" the cov_spectralproj function, "kernel" the cov_kernel function,
"efromovich" the cov_efromovich function and "select" the cov_select function.

* Methods for slm, in "slm-method.R":
The slm function has several associated methods, which are the same as for the lm function.
The available methods are: summary, confint, predict and plot.

* Others functions, in "auxiliary-fun.R":
The package has some auxiliary functions, in particular some predefined kernels for the kernel method of slm function: the
trapeze kernel, the triangle kernel and the rectangular kernel. The user can also define his own kernel and put it in the argument kernel_fonc in the slm function.

* Generative functions, in "generative.R":
The generative_process function generates some stationary processes.
The generative_model function generates some designs.

* Data:
The package contains a dataset "shan". This dataset comes from a study about fine particle pollution in the city of Shanghai. The data are available on the
following website https://archive.ics.uci.edu/ml/datasets/PM2.5+Data+of+Five+Chinese+Cities#.

* References:
E. Caron, J. Dedecker and B. Michel (2019). Linear regression with stationary errors: the R package slm. arXiv preprint arXiv:1906.06583.
https://arxiv.org/abs/1906.06583.
