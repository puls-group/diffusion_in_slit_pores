# Extracting diffusion profiles of particles in slit pores

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10008295.svg)](https://doi.org/10.5281/zenodo.10008295)


This repository contains scripts and cached calculation results for the analysis of anisotropic interface-perpendicular diffusion in film geometries and slit pores.
It most notably contains  a script to calculate the first-order drift correction function `K_B(\gamma)` for bulk-like slabs as presented in the first accompanying publication listed below.
It also contains the scripts to calculate the universal correction functions `R_B(\nu,q)` and `R_LV(\nu,q)` as presented in the second accompanying publication listed below.
The scripts are provided free of charge, as-is to be used in scientific work and educational context. 

If you use the technique presented in the accompanying publication or the scripts in this repository for your own publication, please cite the following publications depending on the aspects of the code you are using.

* For use of the lifetime diffusion analysis method for point-like particles and/or the correction of lifetime based diffusion for drift, please cite:

```
Publication I:
@article{hoellring2023anisotropicI,
title = {Anisotropic molecular diffusion in confinement I: Transport of small particles in potential and density gradients},
journal = {Journal of Colloid and Interface Science},
volume = {650},
pages = {1930-1940},
year = {2023},
issn = {0021-9797},
doi = {https://doi.org/10.1016/j.jcis.2023.07.088},
url = {https://www.sciencedirect.com/science/article/pii/S0021979723013371},
author = {Kevin Höllring and Andreas Baer and Nataša Vučemilović-Alagić and David M. Smith and Ana-Sunčana Smith},
keywords = {Transport coefficient, Diffusion in pores, Porous materials, Anisotropic diffusion, Diffusion at interfaces, Drift and density gradient},
abstract = {Hypothesis
Diffusion in confinement is an important fundamental problem with significant implications for applications of supported liquid phases. However, resolving the spatially dependent diffusion coefficient, parallel and perpendicular to interfaces, has been a standing issue. In the vicinity of interfaces, density fluctuations as a consequence of layering locally impose statistical drift, which impedes the analysis of spatially dependent diffusion coefficients even further. We hypothesise, that we can derive a model to spatially resolve interface-perpendicular diffusion coefficients based on local lifetime statistics with an extension to explicitly account for the effect of local drift using the Smoluchowski equation, that allows us to resolve anisotropic and spatially dependent diffusivity landscapes at interfaces.
Methods and simulations
An analytic relation between local crossing times in system slices and diffusivity as well as an explicit term for calculating drift-induced systematic errors is presented. The method is validated on Molecular Dynamics simulations of bulk water and applied to simulations of water in slit pores.
Findings
After validation on bulk liquids, we clearly demonstrate the anisotropic nature of diffusion coefficients at interfaces. Significant spatial variations in the diffusivities correlate with interface-induced structuring but cannot be solely attributed to the drift induced by local density fluctuations.}
}
```

In this publication, we introduce the concept of lifetime-based diffusion analysis in confinement and deal with the quantification of the effect of drift induced by an effective potential on the observed lifetime statistics.

* For use of the lifetime diffusion analysis of extensive particles, i.e. the universal correction function for either bulk or interface-adjacent slices, or the correction for drift in complex particles,  please cite:

```
Publication II:
@article{hoellring2024anisotropicII,
title = {Anisotropic molecular diffusion in confinement II: A model for structurally complex particles applied to transport in thin ionic liquid films},
journal = {Journal of Colloid and Interface Science},
volume = {657},
pages = {272-289},
year = {2024},
issn = {0021-9797},
doi = {https://doi.org/10.1016/j.jcis.2023.11.137},
url = {https://www.sciencedirect.com/science/article/pii/S0021979723022609},
author = {Kevin Höllring and Andreas Baer and Nataša Vučemilović-Alagić and David M. Smith and Ana-Sunčana Smith},
keywords = {Transport coefficient, Molecular liquids, Diffusion in films, Ionic liquids, Anisotropic diffusion, Diffusion at interfaces},
abstract = {Hypothesis
Diffusion in confinement is an important fundamental problem with significant implications for applications of supported liquid phases. However, resolving the spatially dependent diffusion coefficient, parallel and perpendicular to interfaces, has been a standing issue and for objects of nanometric size, which structurally fluctuate on a similar time scale as they diffuse, no methodology has been established so far. We hypothesise that the complex, coupled dynamics can be captured and analysed by using a model built on the 2-dimensional Smoluchowski equation and systematic coarse-graining.
Methods and simulations
For large, flexible species, a universal approach is offered that does not make any assumptions about the separation of time scales between translation and other degrees of freedom. The method is validated on Molecular Dynamics simulations of bulk systems of a family of ionic liquids with increasing cation sizes where internal degrees of freedom have little to major effects.
Findings
After validation on bulk liquids, where we provide an interpretation of two diffusion constants for each species found experimentally, we clearly demonstrate the anisotropic nature of diffusion coefficients at interfaces. Spatial variations in the diffusivities relate to interface-induced structuring of the ionic liquids. Notably, the length scales in strongly confined ionic liquids vary consistently but differently at the solid–liquid and liquid–vapour interfaces.}
}
```

In this publication, we discuss necessary adaptations of the lifetime based diffusion analysis procedure for the analysis of extensive particles with complex internal degrees of freedom in confinement. 
We extend the drift-correction from the first publication to this more complex scenario and discuss the impact of additional degrees of freedom on lifetime statistics.

## Correcting for drift 

As presented in Publication I, the relative correction for drift induced by an effective potential in confinement can be analytically derived as a function `K_B(\gamma)`, where `\gamma` is a value resulting from a combination of the slope of the effective potential derived from the logarithmic density distribution in the system as well as the chosen slice thickness `L`.
The script to calculate `K_B` can be found at `./scripts/drift_correction_function.py`.


## The universal correction function for bulk-like slabs

### Cached data

In the file `./cache/local_bulk_correction.cache`, we provide tabular pre-calculated values of `R(\nu,q)`.
The values for \nu range from 0.01 to 100 and the range for q is 0.0 to 1.1.

A short explanation of the textual contants of the above mentioned cache file. The most relevant columns are the first three:
* The first column is the value of q=d_mol/L, the relative deformation amplitude 
* The second column is the value of \nu = D_mol / D_\perp, the relative deformation diffusion speed
* The third column is the derived value of R(\nu,q)
* The fourth column is for debugging, denoting how much time has been simulated before applying the long-term extrapolation. It is provided as a relative scale based on the SPM lifetime prediction
* The fifth column denotes how many seconds it took to simulate the EPM system to calculate R(\nu,q)

This description of the file contents can also be found in the file header.

### The script to calculate corrections

We also provide the script written to calculate these values as well as stitch the cached values together for better usability.
It can be found in `./scripts/bulk_correction_function.py`.

#### Dependencies
The script requires the following python modules to work: `scipy`, `numpy`, `matplotlib`
These can be installed using pip. 

#### Correction function
The file contains the function `correction(D_trans: float, slab_width: float, D_displacement: float, displacement: float, skip_cache: bool = False) -> float`, whose parameters we would like to detail. 
It is the function that allows for the stitching of the cache data, but it will also run the adapted Fokker-Planck-Equation (via the function `simulate_fpe_for(rel_displacement: float, rel_D_displacement: float,   dx=0.01, dt=0.1, return_detailed_stats=False)`) and calculate further values of `R(\nu,q)` when no appropriate cache value is found.

The parameters of `correction`:
* `D_trans: float`: The translational diffusion value, denoted by `D_\perp` in the main manuscript
* `slab_width: float`: The slice thickness denoted by `L` in the main manuscript
* `D_displacement: float`: The deformational diffusion value, denoted by `D_mol` in the main manuscript
* `displacement: float`:  The displacement amplitude denoted by `d_mol` in the main manuscript
* `skip_cache: bool = False`: This flag allows to circumvent the cached data and force recalculation of `R(\nu,q)` on each call to the function.

The function `correction` simply calculates the appropriate value of `R(\nu,q)` based on the input parameters and returns the value of `R(\nu,q)`.

#### FPE simulation function

The function `correction` relies on a simulation routine to calculate the result of the FPE presented in the main manuscript to derive `R(\nu,q)`.
This function `simulate_fpe_for(rel_displacement: float, rel_D_displacement: float,   dx=0.01, dt=0.1, return_detailed_stats=False)` has the following parameters:

* `rel_displacement: float`: The relative displacement, denoted by `q` in the main manuscript.
* `rel_D_displacement: float`: The relative deformation diffusion speed, denoted by `\nu` in the main manuscript.   
* `dx=0.01`: The discretisation step in spatial `z` and `s` directions (please refer to the explanation and visualization of the integration domain in the main manuscript)
* `dt=0.1`: The main time step for the RK45 method to simulate the evolution of the density distribution and therefore the statistics of particle lifetimes.
* `return_detailed_stats=False`: A flag to enable more extensive debugging output to be returned.

The function returns a pair of values `(R, t_sim/\tau_SPM)` if `return_detailed_stats` is set to `False`, where `R` denotes the resulting value `R(\nu,q)` and `t_sim/\tau_SPM` denotes the relative scale of the simulated time relative to the lifetime prediction by the SPM. 

If instead `return_detailed_stats` is set to `True`, then a triple of values is returned `(R, t_sim/\tau_SPM, (sim_lt, survival_probability))`, where the first two entries are identical with the case of `return_detailed_stats=False`, but the third entry is a pair of two arrays, one with the time variable `sim_lt` and the second containing the respective survival probabilities `P(sim_lt)` that a particle within the slice at time `t=0` has not left the slice until time `sim_lt`.
Both arrays are of equal length and contain corresponding entries.

## The universal correction function for interface-adjacent slabs

For slices of the system immediatel adjacent to vacuum-like interfaces (i.e. with no strong liquid-interface interactions), we also provide the script to calculate `R_LV(\nu,q)`.
The script can be found at `./scripts/lv_correction_function.py`.