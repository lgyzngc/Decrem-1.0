# Decrem - topologically decoupled metabolic network model
## Installation:
1, install the [**COBRA2.0**](https://opencobra.github.io/cobratoolbox/stable/) toolbaox from palsson lab, and initializing the COBRA tools using command line:
```Bash
initCobraToolbox
```
2, install the matlab toolbox [**SpectraLib_A**](SpectraLib\_A) for clustering analysis;  
3, install the matlab toolbox [**Fast_SNP**](Fast\_SNP) for solving sparse basis vectors;  

## Data preparation:
**Necessary:**<br>
You should first prepare your input files, example is as follows:  
1, the genome-scale metabolic model of given species, e.g. [**bacillus_iYO844**](bacillus/iYO844.mat)  
2, the cofactor metabolites, e.g. [**bacillus_iUO844_cofactors**](bacillus/cofactor.txt)  
3, the elementary enviroment-cell exchange reactions, e.g. [**bacillus_iUO844_exchanges**](bacillus/general_IO_bacillus.txt)  
4, the secreted reactions, e.g. [**bacillus_iUO844_secretion**](bacillus/secrated_bacillus.txt)  
5, the nutrient uptake reactions, e.g. [**bacillus_iUO844_nutrient**](bacillus/nutrient_bacillus.txt)  

**Optional:**<br>
1, experimental 13C measured reaction fluxes for model validation, e.g. [**bacillus_iUO844_13C_flux**](bacillus/intracellularflux_bacillus.txt) 
2, experimental extracellular uptake/secrete rate, e.g. the carbon uptake, O2 rate, CO2 rate, actate or ethonal secretion et.al. which has a same file format with the 13C flux file.
3, the gene knockout list, each row should only store one knockout gene.

## User tutorial:
1, bulid the topologically-decoupled metabolic model  
1, directly running the Decrem_Demo.m for bacillus iYO844 model test or modifing the default input files as user self-defined data.
```Bash
 Decrem_Demo
```
2, bulid the topologically-decoupled metabolic model and predict metabolic fluxes or growth rate using command line:  
```Bash
 [Decrem_solution,Decrem_model,CBM_model] = Decrem(CBM_Model,cofactors,input_nutrient,secretion,general_IO,cluster_num,extraflux(:,[1,i]),intraflux(:,[1,i]),knockout_genes);
``` 
for given CBM_model model,e.g. bacillus iYO844 model or replace all the data in [**bacillus**](bacillus) filefold with user self-defined data  


## pre-trained models:
three reconstructed decoupled metabolic models for Escherichia coli, Saccharomyces cerevisiae and Bacillus subtilis are located in the filefold of 
[**three reconstructed models**](three\sreconstructed\smodels)  

### Reference
**Gaoyang Li et. al. Improved phenotypic predictions of metabolic models by integrating regulatory constraints,BioRxiv,2020**<br>
