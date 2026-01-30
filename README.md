# HII-CHI-Mistry

HII-CHI-Mistry is a collection of python subroutines aimed at the calculation of chemical abundances and physical properties using emission line fluxes from ionised gaseous nebulae. Main website of the code is available [here](https://home.iaa.csic.es/~epm/HII-CHI-mistry.html).

## Requirements

### Python Packages

HII-CHI-Mistry was written for [Python v.2.7](https://www.python.org/download/releases/2.7/), but since its latest versions it is also compatible with [Python v.3](https://www.python.org/download/releases/3.0/). It requires the Python library [NumPy](https://numpy.org).

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

The latest version available is HII-CHI-Mistry [v.6.0](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-opt/HCm_v6.0). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-opt/HCm_v6.0/HCm_v6.0.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Pérez-Montero 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.441.2663P/abstract)).
- Version for Seyferts 2 ([Pérez-Montero et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.2652P/abstract)).
- Version for Low-Ionization Active Galactic Nuclei ([Pérez-Díaz et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.4289P/abstract)).
- Version for the estimation of sulfur abundances in Star-Forming Galaxies ([Pérez-Montero et al. 2025b](https://ui.adsabs.harvard.edu/abs/2025OJAp....851253P/abstract))

### HII-CHI-Mistry ultraviolet range (HCm-uv)

The latest version available is HII-CHI-Mistry-UV [v.5.3](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-uv/HCm-UV_v5.3). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-uv/HCm-UV_v5.3/HCm-UV_v5.3.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Pérez-Montero & Amorín 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.1287P/abstract)).
- Version for Active Galactic Nuclei ([Pérez-Montero et al. 2023](https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.1556P/abstract))

### HII-CHI-Mistry infrared range (HCm-ir)

The latest version available is HII-CHI-Mistry-IR [v.3.2](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-ir/HCm-IR_v3.2). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-ir/HCm-IR_v3.2/HCm-IR_v3.2.readme). The code is described in the following papers:

- Version for Star-Forming Galaxies ([Fernández-Ontiveros et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...652A..23F/abstract)).
- Version for Active Galactic Nuclei ([Pérez-Díaz et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A.115P/abstract)).
- Version for sulfur estimations ([Pérez-Díaz et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...685A.168P/abstract))

### HII-CHI-Mistry effective temperature (HCm-teff)

The latest version available is HII-CHI-Mistry-Tefff [v.5.5](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-teff/HCm-Teff_v5.5). Details on its usage and list of changes can be found [here](https://github.com/Borja-Perez-Diaz/HII-CHI-Mistry/tree/main/HCm-teff/HCm-Teff_v5.5/HCm-Teff_v5.5.readme). The code is described in the following papers:

- Version using plane-parallel or spherical geometry ([Pérez-Montero et al. 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.3322P/abstract)).
- Version using density-bounded models to estimate absorbed photons ([Pérez-Montero et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...643A..80P/abstract)).
- Latest version accepting as input [NII] emission lines, as explained in [Pérez-Montero et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...669A..88P/abstract).


Additionally, a new release of an infrared version of this code has been published, namely HII-CHI-Mistry-Teff-IR [v.3.0](https://github.com/Estallidos/HII-CHI-Mistry/tree/main/HCm-teff/HCm-Teff-IR_v3.0). Details on its usage and list of changes can be found [here](https://github.com/Estallidos/HII-CHI-Mistry/blob/main/HCm-teff/HCm-Teff-IR_v3.0/HCm-Teff-IR_v3.0.readme). The code is described in the following paper:

- Introduction to the code based on infrared emission lines ([Pérez-Montero et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...684A..40P/abstract)).
- Infrared version for the estimation of the AGN components ([Pérez-Montero et al. 2025a](https://ui.adsabs.harvard.edu/abs/2025A%26A...696A.229P/abstract))


## Inputs and outputs


### Optical

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'OII_3727' and 'eOII_3727': emission line ratio [OII] 3726+3729/Hbeta and its error.
- 'NeIII_3868' and 'eNeIII_3868': emission line ratio [NeIII] 3868/Hbeta and its error.
- 'OIII_4363' and 'eOIII_4363': emission line ratio [OIII] 4363/Hbeta and its error.
- 'OIII_4959' and 'eOIII_4959': emission line ratio [OIII] 4959/Hbeta and its error.<sup>[1](#myfootnote1)</sup>
- 'OIII_5007' and 'eOIII_5007': emission line ratio [OIII] 5007/Hbeta and its error.<sup>[1](#myfootnote1)</sup>
- 'NII_5755' and 'eNII_5755': emission line ratio [NII] 5755/Hbeta and its error.
- 'SIII_6312' anf 'eSIII_6312': emission line ratio [SIII] 6312/Hbeta and its error.
- 'NII_6584' and 'eNII_6584': emission line ratio [NII] 6584/Hbeta and its error.
- 'SII_6725' and 'eSII_6725': emission line ratio [SII] 6717+6731/Hbeta and its error. It is possible to give both emission lines separated as 'SII_6717' and 'SII_6731'.
- 'OII_7325' and 'eOII_7325': emission line ratio [OII] 7319+7330/Hbeta and its error.
- 'SIII_9069' and 'eSIII_9069': emission line ratio [SIII] 9069/Hbeta and its error.<sup>[2](#myfootnote2)</sup>
- 'SIII_9532' and 'eSIII_9532': emission line ratio [SIII] 9532/Hbeta and its error.<sup>[2](#myfootnote2)</sup>

<sup><a name="myfootnote1">1</a></sup> It is possible to use only one of the two strong nebular [OIII] emission lines.
<sup><a name="myfootnote2">2</a></sup> It is possible to use only one of the two strong nebular [SIII] emission lines.

**The above emission lines must be reddening corrected**. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- Grid 'i': Index of the grid employe (further details in the .readme file of the code).
- 'O/H' and 'eO/H': estimation of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'S/H' and 'eS/H': estimation of the sulfur abundance 12+log(S/H) and is uncertaintty (only for POPSTAR models).
- 'N/O' and 'eN/O': estimation of the nitrogen-to-oxygen abundance log(N/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Ultraviolet

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'Lya_1216' and 'eLya_1216': emission line Lya HI 1216 and its error.
- 'NV_1239' and 'eNV_1239': emission line NV] 1239 and its error.
- 'CIV_1549' and 'eCIV_1549': emission line CIV 1549 and its error.
- 'HeII_1640' and 'eHeII_1640': emission line HeII 1640 and its error.
- 'OIII_1665' and 'eOIII_1665': emission line OIII] 1665 and its error.
- 'CIII_1909' and 'eCIII_1909': emission line CIII 1909 and its error.
- 'OII_3727' and 'eOII_3727': emission line [OII] 3726+3729 and its error.
- 'NeIII_3868' and 'eNeIII_3868': emission line [NeIII] 3868 and its error.
- 'OIII_4363' and 'eOIII_4363': emission line [OIII] 4363and its error.
- 'Hb_4861' and 'eHb_4861': emission line Hb HI 4861 and its error.
- 'OIII_5007' and 'eOIII_5007': emission line [OIII] 5007 and its error.

**The above emission lines must be reddening corrected and in the same units**. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-uv-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- Grid 'i': Index of the grid employe (further details in the .readme file of the code).
- 'O/H and 'eO/H': estimation of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'C/O and 'eC/O': estimation of the carbon-to-oxygen abundance log(C/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Infrared

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- 'HI_4m' and 'eHI_4m': emission line HI 4.05 mic and its error.
- 'ArII_7m' and 'eArII_7m': emission line [ArII] 6.98 mic and its error.
- 'HI_7m' and 'eHI_7m': emission line HI 7.46 mic and its error.
- 'ArV_8m' and 'eArV_8m': emission line [ArV] 7.90 mic and its error.
- 'ArIII_9m' and 'eArIII_9m': emission line [ArIII] 8.99 mic and its error.
- 'SIV_10m' and 'eSIV_10m': emission line [SIV] 10.5 mic and its error.
- 'HI_12m' and 'eHI_12m': emission line HI 12.4 mic and its error.
- 'NeII_12m' and 'eNII_12m': emission line [NeII] 12.8 mic and its error.
- 'ArV_13m' and 'eArV_13m': emission line [ArV] 13.1 mic and its error.
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
- 'S/H' and 'eS/H': estimation of the sulfur abundance 12+log(S/H) and its uncertainty.
- 'N/O' and 'eN/O': estimation of the nitrogen-to-oxygen abundance log(N/O) and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

### Effective temperature (optical version)

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- '12logOH' and 'e12logOH': if known, oxygen abundance 12+log(O/H) and its error.<sup>[3](#myfootnote3)</sup>
- 'OII_3727' and 'eOII_3727': emission line [OII] 3727 and its error.
- 'OIII_4959' and 'eOIII_3727': emission line [OIII] 4959 and its error.<sup>[4](#myfootnote4)</sup>
- 'OIII_5007' and 'eOIII_5007': emission line [OIII] 5007 and its error.<sup>[4](#myfootnote4)</sup>
- 'SII_6725' and 'eSII_6725': emission line [SII] 6717+6731 and its error. It is possible to give both emission lines separated as 'SII_6717' and 'SII_6731'.
- 'SIII_9069' and 'eSIII_9069': emission line [SIII] 9069 and its error.<sup>[5](#myfootnote5)</sup>
- 'SIII_9532' and 'eSIII_9532': emission line [SIII] 9532 and its error.<sup>[5](#myfootnote5)</sup>
- 'HeI_4471' and 'eHeI_4471': emission line HeI 4471 and its error.
- 'HeI_5876' and 'eHeI_5876': emission line HeI 5876 and its error.
- 'HeII_4686' and 'eHeII_4686': emission line HeII 4686 and its error.
- 'ArIV_4730' and 'eAIV_4740': emission line [ArIV] 4740 and its error.
- 'ArIII_7135' and 'eAIII_7135': emission line [ArIII] 7135 and its error.
- 'NII_6584' and 'eNII_6584': emission line [NII] 6584 and its error.

<sup><a name="myfootnote3">3</a></sup> If oxygen abundances are unknown, HCm will estimate them.
<sup><a name="myfootnote4">4</a></sup> It is possible to use only one of the two strong nebular [OIII] emission lines.
<sup><a name="myfootnote5">5</a></sup> It is possible to use only one of the two strong nebular [SIII] emission lines.

The above emission lines must be reddening corrected. In case no errors are provided, it is advisable not using MonteCarlo iterations at all. All comments must be placed at the beginning of the document preceded by the symbol "#".

The output file is also a text format, named with the original name of the input file and followed by the extension "_hcm-teff-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- 'O/H' and 'eO/H': estimation or input value of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'Teff' and 'eTeff': estimation of the effective temperature in K and its uncertainty.<sup>[6](#myfootnote6)</sup>
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

<sup><a name="myfootnote6">6</a></sup> If required, columns 'Teff' and 'eTeff' are replaced by 'f_abs' and 'ef_abs' (fraction of absorbed photons).

### Effective temperature (infrared version)

The input file must be written in text format with a first row of labels indicating some or all of the following columns:

- 'ID': identification for each row.
- '12logOH' and 'e12logOH': if known, oxygen abundance 12+log(O/H) and its error.<sup>[7](#myfootnote7)</sup>
- 'ArII_7m' and 'eArII_7m': [ArII] 6.98 mic and its error
- 'ArV_8m' and 'eArV_8m': [ArV] 7.90 mic and its error
- 'ArIII_9m' and 'eArIII_9m': [ArIII] 8.99 mic and its error
- 'SIV_10m' and 'eSIV_10m': [SIV] 10.5 mic and its error
- 'NeII_12m' and 'eNeII_12m': [NeII] 12.8 mic and its error
- 'ArV_13m' and 'eArV_13m': [ArV] 13.1 mic  and its error
- 'NeV_14m' and 'eNeV_14m': [NeV] 14.3 mic  and its error
- 'NeIII_15m' and 'eNeIII_15m': [NeIII] 15.5 mic  and its error
- 'SIII_18m' and 'eSIII_18m': [SIII] 18.7 mic and its error
- 'NeV_24m' and 'eNeV_24m': [NeV] 24.2 mic  and its error
- 'OIV_25m' and 'eOIV_25m': [OIV] 25.9 mic  and its error
- 'SIII_33m' and 'eSIII_33m': [SIII] 33.7 mic and ist error
- 'OIII_52m' and 'eOIII_52m': [OIII] 52 mic and its error
- 'NIII_57m' and 'eNIII_57m': [NII] 57 mic and its error
- 'OIII_88m' and 'eOIII_88m': [OIII] 88 mic and its error
- 'NII_122m' and 'eNII_122m': [NII] 122 mic and its error
- 'NII_205m' and 'eNII_205m': [NII] 205 mic and its error

<sup><a name="myfootnote7">7</a></sup> If oxygen abundances are unknown, HCm will estimate them.

The above emission lines must be reddening corrected. If no information exists about a certain column it must typed as zero or not introduced in the input file. If the error is not known or if it is not going to be taken into account in the calculations,  it is advisable not using MonteCarlo iterations at all. Regarding lines the routine will only provide a  solution if at least one low-to-high excitation is given (e.g. [NeII] and [NeIII] and/or [SIII] and [SIV]). If only two low-excitation or high-excitation lines  are given the program will provide 0 values in the results.

The output file is also a text format, named with the original name of the input file and followed by the extension "_HCm-Teff-IR-output.dat". The first columns will show the information provided in the input file. In case ID column is missing, the code will automatically assign a cardinal to each row. The last seven columns show:

- 'O/H' and 'eO/H': estimation or input value of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'Teff' and 'eTeff': estimation of the effective temperature in K and its uncertainty.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

If the AGN version is chosen instead, the outputs are:

- 'O/H' and 'eO/H': estimation or input value of the oxygen abundance 12+log(O/H) and its uncertainty.
- 'a_OX' and 'ea_OX': estimation of the a(OX) slope in the AGN SED.
- 'U' and 'eU': estimation of the ionisation parameter log(U) and its uncertainty.

## Non-interactive mode

The latest versions of the code allow the user to implement a non-interactive mode that allows the code to be implemented within any other python script. All the details in the selection of the model must be modified within the python script. For example, in the optical version the variables that must be modified are:

```python
#################################################
###### NON INTERACTIVE OR INTERACTIVE CODE ######
#################################################

interactive = True #Change this value to False to run the code in non-interactive mode

#Questions (inputs fromt terminal)
question1 = interactive #Question to select the grids of models
question2 = interactive #Question for the value of alphaOX (only required if "sed" is set to 3)
question3 = interactive #Question for the value of efrac (only required if "sed" is set to 3)
question4 = interactive #Question required to introduce a particular file by the user (only required if "sed" is set 4)
question6 = interactive #Question to use or not interpolation for the grids
question7 = interactive #Question regarding the constrain law to be assumed by the code
question8 = interactive #Question required to tntroduced a particular file by the user (only if "const" is set to 6)

#Set values to operate
if question1 == False:
   sed = 1 #Choose grid of model: (1) for POPSTAR models; (2) for BPASS models; (3) for AGN models; (4) models introduced by the user (will require file). Replace value for the given option
if question2 == False and sed == 3:
   alpha = 2 #Choose value of alpha_OX: (1) for alpha_OX = -0.8; (2) for alpha_OX = -1.2. Replace value with any of these options
if question3 == False and sed == 3:
   efrac = 1 #Choose value for the stopping criteria, i.e., the fraction of free electrons: (1) for efrac = 0.02; (2) for efrac = 0.98
if question4 == False and sed == 4:
   new_library = 'Name_of_the_file' #Introduced name of the file with the grids of models. It must be located under the folder "Libraries_ir"
if question6 == False:
   inter = 0 #Choose value to perform interpolation: (0) no interpolation; (1) interpolation.  Replace value for the given option
if question7 == False:
   const = 1 #Choose value for the constraint laws between O/H, N/O and U that must be assumed when the code does not have enough information: (1) constraints obtained for Star-Forming Galaxies; (2) constraints obtained for Extreme Emission Line Galaxies; (3) constraints between N/O and O/H obtained for AGNs, without restriction in the ionization parameter; (4) constraints between N/O and O/H obtained for AGNs, and log(U) > -2.5; (5) constraints between N/O and O/H obtained for AGNs, and log(U) < -2.5; (6) constraint law introduced by the user (will required a file).  Replace value for the given option
if question8 == False and const == 6:
   new_const = 'Name_of_the_file' #Introduced name of the file with the constraint laws. It must be located under the folder "Constraint"
```

A similar structure is presented in the rest of flavours of the code.

## Support

This program has been made thanks to the financial support from the Spanish AYA project Estallidos.

## Contact

Further questions, comments and suggestions are welcome to:

- Enrique Pérez-Montero | Pronouns: he/his | Affiliation: [IAA-CSIC](https://www.iaa.csic.es) | Publications: [List](https://ui.adsabs.harvard.edu/search/q=%20%20author%3A%22Perez-Montero%2C%20E.%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) | Mail: [epm[at]iaa.es](mailto:epm@iaa.es)
- Borja Pérez-Díaz | Pronouns: he/his | Affiliation: [OAR-INAF](https://www.oa-roma.inaf.it) | Publications: [List](https://ui.adsabs.harvard.edu/search/q=%20%20author%3A%22Perez-Diaz%2C%20Borja%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) | Mail: [borja.perezdiaz[at]inaf.it](mailto:borja.perezdiaz@inaf.it)
- Juan Antonio Fernández-Ontiveros | Pronouns: he/his | Affiliation: [CEFCA](https://www.cefca.es) | Publications: [List](https://ui.adsabs.harvard.edu/search/q=%20%20author%3A%22Fernandez-Ontiveros%22&sort=date%20desc%2C%20bibcode%20desc&p_=0) | Mail: [j.a.fernandez.ontiveros[at]gmail.com](mailto:j.a.fernandez.ontiveros@gmail.com)
