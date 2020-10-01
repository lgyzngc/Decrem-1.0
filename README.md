# Decrem - topologically decoupled metabolic network model
## Installation:
1, install the [**COBRA2.0**](https://opencobra.github.io/cobratoolbox/stable/) toolbaox from palsson lab.
2, initializing the COBRA tools using command line:
```Bash
initCobraToolbox
```
3, add java package of [**jartest**](jartest) into matlab for simple directed cycle identification using the command:   
```Bash
 javaclasspath
.\jartest\lib\jgraph-5.13.0.0.jar                                                   
.\jartest\lib\jgrapht-core-0.9.0.jar                                                
.\jartest\lib\jgrapht-demo-0.9.0.jar                                                
.\jartest\lib\jgrapht-ext-0.9.0.jar                                                 
.\jartest\lib\jgrapht-ext-0.9.0-uber.jar                                            
.\jartest\lib\jgraphx-2.0.0.1.jar  
.\jartest\simplecyclesofklength.jar 
```
4, install the matlab toolbox [**SpectraLib_A**](SpectraLib\_A) for clustering analysis;  
5, install the matlab toolbox [**Fast_SNP**](Fast\_SNP) for solving sparse basis vectors;  

## Data preparation:
**Necessary:**<br>
Your should first prepare your input files, example is as follows:  
1, the genome-scale metabolic model of given species, i.e. [**bacillus_iYO844**](bacillus/iYO844.mat)  
2, the elementary enviroment-cell exchange reactions, i.e. [**bacillus_iUO844_exchanges**](bacillus/general_IO_bacillus.txt)  
3, the nutrient uptake reactions, i.e. [**bacillus_iUO844_nutrient**](bacillus/nutrient_bacillus.txt)  

**Optional:**<br>
1, experimental 13C measured reaction fluxes for model validation, i.e. [**bacillus_iUO844_13C_flux**](bacillus/intracellularflux_bacillus.txt)  
2, simple cycle-derived reaction similarity metrix, i.e. [**bacillus_iUO844_similarity**](bacillus/similarity_matrix_5len_rec4.txt)  


## User tutorial:
1, bulid the topologically-decoupled metabolic model  
```Bash
 test_reconstruction
```
for bacillus iYO844 model, or using command line:  
```Bash
 reconstructed_model = decoupledModelConstruct(model,cofactor_path,secrated_path,nutrient_path,general_IO_path)
 %%model: metabolic models with mat format
 %%cofactor_path: the file path of cofactors
 %%secrated_path: the file path of secrated reactions
 %%nutrient_path: the file path of nutrient reactions
 %%general_IO_path: the file path of general IO reactions
```
for usr self-defined data

2, predict metabolic fluxes or growth rate using command line:  
```Bash
 Decrem
``` 
for bacillus iYO844 model, or replace all the data in bacillus filefold with usr self-defined data  

### Reference
**Gaoyang Li et. al. Improved phenotypic predictions of metabolic models by integrating regulatory constraints,BioRxiv,2020**<br>
