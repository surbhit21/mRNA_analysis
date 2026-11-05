
This is the read me for the repository containing the scripts used to analyse mRNA localization and generate figures for the paper
"An integrative model of AMPA receptor trafficking reveals the central contribution of local translation in subtype-specific kinetics"
Surbhit Wagle, Maximilian K. Kracht, Anne Bührke, Amparo Acker-Palmer, Nataliya Kraynyukova, Anne-Sophie Hafner,  Erin M. Schuman,  Tatjana Tchumatchenko
doi: https://doi.org/10.1101/2025.02.08.637220

# 1.System requirements 
The requirement for running this code includes: Python and jupyter notebook. The code was developed on Apple MacBook Pro M1 running Sonoma 14.5 os with a RAM of 32 GB and
storage of 1 TB. As some of the figure requires downloading simulation data, please insure you have enough RAM and storage space.

# 2.Pre-requisite

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install following packages.

```bash
asteval==0.9.28
contourpy==1.0.6
cycler==0.11.0
fonttools==4.38.0
future==0.18.2
kiwisolver==1.4.4
lmfit==1.0.3
matplotlib==3.6.2
numpy==1.23.4
packaging==21.3
pandas==1.5.1
patsy==0.5.3
Pillow==9.3.0
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.6
scikit-posthocs==0.7.0
scipy==1.9.3
seaborn==0.12.1
six==1.16.0
statsmodels==0.13.5
uncertainties==3.1.7
```
To install the above-listed packages, run: 
```bash
pip install -r requirement.txt 
```
Or you can choose to install them separately
# 3.Usage

Download or clone this repository and install the dependencies. The table below identifies the analysis entry points for the figure panels. Plotting source data are provided in `Figure_data`; rerunning the original analyses additionally requires the corresponding raw or intermediate input files which are not necessarily included in the repository (see Data Sources).

Run Python commands from the repository root. For notebooks, start Jupyter in `Scripts` so that relative paths such as `./Protein_data` resolve correctly:

```bash
cd Scripts
jupyter notebook
```

Open the relevant notebook and run its cells in order with the required inputs available.

| Figure panels | Analysis script or notebook | Analysis |
| --- | --- | --- |
| 1B–C | `python Scripts/mRNA_Analysis_v3.py -m Gria2` | Somatic/dendritic fractions and dendritic distribution of *Gria2* versus *Camk2a* mRNA |
| 1E | `python Scripts/CNIH2_protein_analysis.py` | MAP2-normalized total GluA2 distribution |
| 1F–H | `python Scripts/DifferentParams.py` | Model parameter comparisons |
| 2B–C | `python Scripts/HaloGluA2_analysis.py` | Surface/intracellular GluA2 ratios for the two labeling schemes |
| 2F | `python Scripts/SE_HaloTag.py` | Spine-to-shaft GluA2 enrichment |
| 4B–C | `python Scripts/mRNA_Analysis_v3.py -m CNIH2` | Somatic/dendritic fractions and dendritic distribution of *Cnih2* versus *Camk2a* mRNA |
| 4E–F | `Scripts/CNIH2_LTP.ipynb` | Nascent CNIH2 in basal and cLTP conditions, in soma and dendrites |
| 4H–I | `Scripts/CNIH2_GluA2.ipynb` | Nascent surface GluA1 and GluA2 under control and *Cnih2*-targeting shRNA |
| S1B–C | `python Scripts/mRNA_Analysis_v3.py -m Gria1` | Somatic/dendritic fractions and dendritic distribution of *Gria1* versus *Camk2a* mRNA |
| S2 | `python Scripts/DifferentParams.py` | Model parameter comparisons |
| S3 | `python Scripts/CNIH2_protein_analysis.py` | CNIH2/GluA2 protein analysis |
| S4A–B | `python Scripts/HaloGluA2_analysis.py` | Normalized surface and intracellular GluA2 fluorescence profiles, 100 µm selection |
| S5A–H | `python Scripts/HaloGluA2_analysis.py` | Binned surface/intracellular ratios for 75, 100, 125 and 150 µm selections, respectively |
| S6C | `python Scripts/Spine_Enrichment_comparison.py` | Antibody versus HaloTag surface GluA2 synaptic enrichment |
| S10B–D | `Scripts/CNIH2_new_synthesis.ipynb` | Nascent CNIH2 in Met versus AHA conditions and AHA soma/dendrite fractions |
| S10F | `python Scripts/CNIH2_protein_analysis.py` | MAP2-normalized CNIH2 distribution |
| S11 | `Scripts/CNIH2_validation.ipynb` | Validation of *Cnih2* knockdown |
| S12B | `Scripts/CNIH2_GluA2_cLTP.ipynb` | Control shRNA, *Cnih2*-targeting shRNA and *Cnih2*-targeting shRNA with cLTP |
| S15B–C | `python Scripts/SE_HaloTag.py` | Spine density before and after filtering, per cell and per dendrite |


## Figure 5 simulation plots

Download the simulation data from [CNIH2_simulations](https://gin.g-node.org/surbhitw/CNIH2_simulations). Set `input_folder` in `Scripts/LoadNPlot.py` to the downloaded data location. The dataset identifiers previously used for the three panels are:

```bash
# Figure 5A
python Scripts/LoadNPlot.py -d 04_19_2026_18_10_07

# Figure 5B
python Scripts/LoadNPlot.py -d 04_19_2026_18_57_10

# Figure 5C
python Scripts/LoadNPlot.py -d  04_19_2026_20_55_28
```

# 4.Contribution

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home)

# 4.License

[MIT](https://choosealicense.com/licenses/mit/)


# 5.Data Sources
All data used for plotting are provided in the `Figure_data` folder. Raw data required to rerun some of the analyses can be made available upon reasonable request.

# 6.Acknowledgment
We would like to thank all Tchumathenko, Hafner, Acker-Palmer, Schuman and Choquet lab members for fruitful discussions and helpful feedback.

This study was supported by the University of Bonn Medical Center (SW, CB, NK, TT), the University of Mainz Medical Center (SW, TT), the German Research Foundation via CRC1080 (SW, MKK, AAP, EMS, TT), CRC1507 (AAP), the Donders Institute for Brain, Cognition and Behaviour and Faculty of Science, Radboud University Nijmegen Netherlands (AH). This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (‘MolDynForSyn’, grant agreement No. 945700) (TT) & (‘MemCode’, grant agreement No. 101076961) (AH). AH also received support from the EMBO long-term postdoctoral fellowship (ALTF 1095-2015) and the Alexander von Humboldt Foundation (FRA-1184902-HFST-P).
