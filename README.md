# INLAspacetime

This is a R package to implement spatial and spatio-temporal models
taking use to the `cgeneric` interface in the INLA package. 
This interface is a way to implement models by writing `C`
code to build the precision matrix compiling it so that 
INLA can use it internally.

# So far we have implemented 

1. some of the models presented in 
https://arxiv.org/abs/2006.04917

2. the barrier model proposed in 
https://doi.org/10.1016/j.spasta.2019.01.002

# Install with 
```
remotes::install_github("eliaskrainski/INLAspacetime")
```

# See the vignettes for examples
