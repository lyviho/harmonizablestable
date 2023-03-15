## Non-ergodic statistics and spectral density estimation for stationary real harmonizable symmetric alpha-stable moving averages
### Ly Viet Hoang and Evgeny Spodarev
Institute of Stochastics, Ulm University, Ulm Germany

https://doi.org/10.48550/arXiv.2209.04315

### Multiple independent paths
Matlab files for the simulation and estimation with L independent paths
can be found in the folder `mult_paths`. Main file is `sd_est_alphasine.m`,
the rest are auxiliary functions needed for the path simulation of the process, 
estimation of the codifference function as well as the alpha-sine inversion. 

For theoretical details see Sections 2.2. and 5.1.

### Single path
Matlab and R files in `single_path`. The matlab skript `per_freq_est.m` generates a 
path of the process X and estimates the underlying i.i.d. frequencies Z_k with spectral
density as their pdf. The results are saved as ´filename.mat´.
`KDE.Rmd` imports this `.mat` file and computes the kernel density estimator. 
We use R because of its already integrated density estimation with various 
kernel functions and bandwidth methods.

Examples are included in the subfolder `matdata` for examples f1 and f2 in 
`per_freq_est.m`, e.g. `f1a150_2.mat` corresponds to the frequency estimates for
example f1 with alpha=1.5 and sample size 10^4. See `KDE.Rmd` for more details. 

Details in Sections 4 and 5.2. of the paper. 
