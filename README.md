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


## Inputs and outputs


### Optical

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'OII_3727' and 'eOII_3727': emission line ratio [OII] 3727/Hbeta and its error.
- 'NeIII_3868' and 'eNeIII_3868': emission line ratio [NeIII] 3868/Hbeta and its error.
- 'OIII_4363' and 'eOIII_4363': emission line ratio [OIII] 4363/Hbeta and its error.
- 'OIII_4959' and 'eOIII_4959': emission line ratio [OIII] 4959/Hbeta and its error{1}.
- 'OIII_5007' and 'eOIII_5007': emission line ratio [OIII] 5007/Hbeta and its error{1}.
- 'NII_6584' and 'eNII_6584': emission line ratio [NII] 6584/Hbeta and its error.
- 'SII_6725' and 'eSII_6725': emission line ratio [SII] 6717+6731/Hbeta and its error. It is possible to give both emission lines separated as 'SII_6717' and 'SII_6731'.

{1} **Note**: It is possible to use only one of the two strong nebular [OIII] emission lines.

The above emission lines must be reddening corrected. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- Grid 'i': Index of the grid employe (further details in the .readme file of the code).
- 'O/H' and 'eO/H': estimation of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'N/O' and 'eN/O': estimation of the nitrogen-to-oxygen abundance log(N/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Ultraviolet

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'Lya_1216' and 'eLya_1216': emission line Lay HI 1216 and its error.
- 'CIV_1549' and 'eCIV_1549': emission line CIV 1549 and its error.
- 'HeII_1640' and 'eHeII_1640': emission line HeII 1640 and its error.
- 'OIII_1665' and 'eOIII_1665': emission line OIII] 1665 and its error.
- 'CIII_1909' and 'eCIII_1909': emission line CIII 1909 and its error.
- 'Hb_4861' and 'eHb_4861': emission line Hb HI 4861 and its error.
- 'OIII_5007' and 'eOIII_5007': emission line [OIII] 5007 and its error.

The above emission lines must be reddening corrected and in the same units. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-uv-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- Grid 'i': Index of the grid employe (further details in the .readme file of the code).
- 'O/H and 'eO/H': estimation of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'C/O and 'eC/O': estimation of the carbon-to-oxygen abundance log(C/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Infrared

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'HI_4m' and 'eHI_4m': emission line HI 4.05 mic and its error.
- 'HI_7m' and 'eHI_7m': emission line HI 7.46 mic and its error.
- 'SIV_10m' and 'eSIV_10m': emission line [SIV] 10.5 mic and its error.
- 'HI_12m' and 'eHI_12m': emission line HI 12.4 mic and its error.
- 'NeII_12m' and 'eNII_12m': emission line [NeII] 12.8 mic and its error.
- 'NeV_14m' and 'eNeV_14m': emission line [NeV] 14.3 mic and its error.
- 'NeIII_15m' and 'eNeIII_15m': emission line [NeIII] 15.5 mic and its error.
- 'SIII_18m' and 'eSIII_18m': emission line [SIII] 18.7 mic and its error.
- 'NeV_24m' and 'eNeV_24m': emission line [NeV] 24.2 mic and its error.
- 'OIV_26m' and 'eOIV_26m': emission line [OIV] 25.9 mic and its error.
- 'SIII_33m' and 'eSIII_33m': emission line [SIII] 33.7 mic and its error.
- 'OIII_52m' and 'eOIII_52m': emission line [OIII] 52 mic and its error.
- 'NIII_57m' and 'eNIII_57m': emission line [NIII] 57 mic and its error.
- 'OIII_88m' and 'eOIII_88m': emission line [OIII] 88 mic and its error.
- 'NII_122m' and 'eNII_122m': emission line [NII] 122 mic and its error.
- 'NII_205m' and 'eNII_205m': emission line [NII] 205 mic and its error.

The above emission lines must be given in the same units. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-ir-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- Grid 'i': Index of the grid employe (further details in the .readme file of the code).
- 'O/H' and 'eO/H': estimation of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'N/O' and 'eN/O': estimation of the nitrogen-to-oxygen abundance log(N/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Effective temperature

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- '12logOH' and 'e12logOH': if known, oxygen abundance 12+log(O/H) and its error{1}.
- 'OII_3727' and 'eOII_3727': emission line ratio [OII] 3727/Hbeta and its error.
- 'OIII_4959' and 'eOIII_3727': emission line ratio [OIII] 4959/Hbeta and its error{2}.
- 'OIII_5007' and 'eOIII_5007': emission line ratio [OIII] 5007Hbeta and its error{2}.
- 'SII_6725' and 'eSII_6725': emission line ratio [SII] 6717+6731/Hbeta and its error. It is possible to give both emission lines separated as 'SII_6717' and 'SII_6731'.
- 'SIII_9069' and 'eSIII_9069': emission line ratio [SIII] 9069/Hbeta and its error{3}.
- 'SIII_9532' and 'eSIII_9532': emission line ratio [SIII] 9532/Hbeta and its error{3}.
- 'HeI_4471' and 'eHeI_4471': emission line ratio HeI 4471/Hbeta and its error.
- 'HeI_5876' and 'eHeI_5876': emission line ratio HeI 5876/Hbeta and its error.
- 'HeII_4686' and 'eHeII_4686': emission line ratio HeII 4686/Hbeta and its error.

{1} **Note**: If oxygen abundances are unknown, Him will estimate them.
{2} **Note**: It is possible to use only one of the two strong nebular [OIII] emission lines.
{3} **Note**: It is possible to use only one of the two strong nebular [SIII] emission lines.

The above emission lines must be reddening corrected. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-ir-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- 'O/H' and 'eO/H': estimation or input value of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'Teff' and 'eTeff': estimation of the effective temperature in K and its uncertainty{4}.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

{4} **Note**: If required, columns 'Teff' and 'eTeff' are replaced by 'f_abs' and 'ef_abs' (fraction of absorbed photons).

## Contact

Further questions, comments and suggestions are welcome to:

- Enrique Pérez-Montero | Pronouns: he/his | Affiliation: [IAA-CSIC](https://www.iaa.csic.es) | Publications: [List](https://ui.adsabs.harvard.edu/search/q=%20%20author%3A%22Perez-Montero%2C%20E.%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) | Mail: [epm[at]iaa.es](mailto:epm@iaa.es)
- Borja Pérez-Díaz | Pronouns: he/his | Affiliation: [IAA-CSIC](https://www.iaa.csic.es) | Publications: [List](https://ui.adsabs.harvard.edu/search/q=%20%20author%3A%22Perez-Diaz%2C%20Borja%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) | Mail: [bperez[at]iaa.es](mailto:bperez@iaa.es)
