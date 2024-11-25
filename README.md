# AMPA receptor mRNA figures


# Project

This project contains the scripts used to analyse mRNA localization and generate figures for the paper Wagle et al. 

## Pre-requisite

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
pip install -r requirment.txt 
```
Or you can choose to install them separately
## Usage
Please note that for legacy reasons, we have devided the code base into two repositories.
Here is the link to the second repository. 
To reproduce the figures from Wagle et. al., please do the following
(in case you do not see the figure listed below, please check the second repository)
For Fig 2 B and C
```python
python Scripts/mRNA_Analysis_v3.py -m Gria2
```

For Fig 2E
```python
python Scripts/CNIH2_protein_analysis.py
```

For Fig 2 F-H
```python
python Scripts/DifferentParams.py
```

For Fig 5 B and C

```python
python Scripts/mRNA_Analysis_v3.py -m CNIH2
```

For Fig 5 E and F

```python
CNIH2_LTP.ipynb
```

For Fig 5 H and I

```python
Run the CNIH2_LTP.ipynb
```

For Fig S1 B-C
```python
python Scripts/mRNA_Analysis_v3.py -m Gria1
```

For S2
```python
python Scripts/DifferentParams.py
```

For Fig S3

```python
python Scripts/CNIH2_protein_analysis.py
```

Fig S10 B-D

```python
Run the CNIH_new_synthesis.ipynb
```

For Fig S10F

```python
python Scripts/CNIH2_protein_analysis.py
```

For Fig S11
```python
Run the CNIH_validation.ipynb
```

For Fig S12
```python
python Script/PlotMedianCNIH2.py
```

## Contributing

This code was developed by [Surbhit Wagle](https://sites.google.com/view/surbhitwagle/home)

## License

[MIT](https://choosealicense.com/licenses/mit/)
## Support
This work is supported by the [CRC 1080](https://www.crc1080.com) funded by DFG .


## Authors and acknowledgment
We would like to thank all our lab members for fruitful discussions and helpful feedback. This study was supported by the University of Bonn Medical Center (SW, NK, TT), University of Mainz medical center (SW, TT), the German Research Foundation via CRC1080 (SW, TT), the Donders Institute for Brain, Cognition and Behaviour and Faculty of Science, Radboud University Nijmegen Netherlands (AH). This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (‘MolDynForSyn’, grant agreement No. 945700) (TT) & (‘MemCode’, grant agreement No. 101076961) (AH). AH also received support by the EMBO long-term postdoctoral fellowship (ALTF 1095-2015) and the Alexander von Humboldt Foundation (FRA-1184902HFST-P).

## License
For open source projects, say how it is licensed.

## Project status
We are curretly still working on updating the scripts for small changes in the plots and any new analysis that might come up.
