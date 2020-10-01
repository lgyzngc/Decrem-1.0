# Decrem - topologically-decoupled metabolic network
## Installation:
1, install the [**COBRA2.0**](https://opencobra.github.io/cobratoolbox/stable/) toolbaox from palsson lab.
2, initializing the COBRA tools using command line:
```Bash
initCobraToolbox
```
2, add java package of jartest into matlab using the command:   
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

3, setup the matlab package 'SpectraLib_A' for clustering analysis;  
4, setup the matlab package 'Fast_SNP' for solving sparse basis vectors;  
5, test the 'Decrem.m';  
