# HII-CHI-Mistry

HII-CHI-Mistry is a collection of python subroutines aimed at the calculation of chemical abundances and physical properties using emission line fluxes from ionised gaseous nebulae. 

## Requirements

### Python Packages

HII-CHI-Mistry was written for [Python v.2.7](https://www.python.org/download/releases/2.7/), but since its latest versions it is also compatible with [Python v.3.](https://www.python.org/download/releases/3.0/). It requires the Python library [NumPy](https://numpy.org).

### Libraries

HII-CHI-Mistry uses grid of models (libraries) to perform the estimations. Each version of code has its own libraries, created from photoionization models computed with [Cloudy v.17](https://gitlab.nublado.org/cloudy/cloudy).

Most recent versions of the codes allow the user to introduce a different file. However, this file introduced must be located under the corresponding folder and with the proper format (see documentation for each code).

In any case, **defaults libraries must not be removed**.


## Versions

HII-CHI-Mistry presents different variations depending on the spectral range analysed and the analysis to be performed with the input data. These are the available options:

- HII-CHI-Mistry for optical emission lines ([HCm](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-opt)).
- HII-CHI-Mistry for ultraviolet emission lines ([HCm-uv](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-uv)).
- HII-CHI-Mistry for infrared emission lines ([HCm-ir](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-ir)).
- HII-CHI-Mistry for equivalent effective temperature ([HCm-teff](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-teff)).

### HII-CHI-Mistry optical range (HCm-opt)

The latest version available is HII-CHI-Mistry [v.5.2](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-opt/HCm_v5.2). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-opt/HCm_v5.2/HCm_v5.2.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Pérez-Montero 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.441.2663P/abstract)).
- Version for Seyferts 2 ([Pérez-Montero et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.2652P/abstract)).
- Version for Low-Ionization AGN ([Pérez-Díaz et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.4289P/abstract)).

### HII-CHI-Mistry ultraviolet range (HCm-uv)

The latest version available is HII-CHI-Mistry-UV [v.4.2](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-uv/HCm-UV_v4.2). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-uv/HCm-UV_v4.2/HCm-UV_v4.2.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Pérez-Montero & Amorín 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1287P/abstract)).

### HII-CHI-Mistry infrared range (HCm-ir)

The latest version available is HII-CHI-Mistry-IR [v.2.2](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-ir/HCm-IR_v2.2). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-ir/HCm-IR_v2.2/HCm-IR_v2.2.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Fernández-Ontiveros et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...652A..23F/abstract)).

### HII-CHI-Mistry effective temperature (HCm-teff)

The latest version available is HII-CHI-Mistry-Tefff [v.5.1](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-teff/HCm-Teff_v5.1). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-teff/HCm-Teff_v5.1/HCm-Teff_v5.1.readme). The code is described in the following papers:

- Version using plane-parallel or spherical geometry ([Pérez-Montero et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.3322P/abstract)).
- Version using density-bounded models to estimate absorbed photons ([Pérez-Montero et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...643A..80P/abstract))














