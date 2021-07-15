# 2021_cellsystems

This repository stores the code used to analyze data and generate plots for "Assembling stable syntrophic Escherichia coli communities by comprehensively identifying beneficiaries of secreted goods"

1. strictAndPermissiveAuxotrophy.m 
  analyzes data of deletion strains individually grown in nutrient poor media or nutrient poor media conditioned by a wild-type strain. It calculates a collaboration score based on the growth of an ideal auxotroph. It also compares the growth of individual strains to strains grown as a pooled culture
  it needs to be run with:
    cmAUC.mat
      holds the names of the deletion strains used and the area under the growth curve values for deletion strains grown in nutrient poor media or in nutrient poor media conditioned by a wild-type strain   
   supplementary table 1, sheet "Cond. media 5X10^9 cells per mL" 
      holds the fold change values of pooled deletion strains grown in nutrient poor media or in nutrient poor media conditioned by a wild-type strain culture at a density of 5X10^9 cells per mL

2. pairAnalysis.m 
  analyzes absorbance data from co-cultures to identify pairs that improve their growth in co-culture relative to individual cultures. It uses the area under the growth curve to assess cooperation
  needs to be run with:
    complementation.mat
      holds the absorbance measurements of 1444 auxotroph pairs grown in nutrient poor media 
      
3. groupInteractionNetwork.m 
  takes the pairwise interactions, pools them together by their functional group and builds a network. It makes an xlsx file that can be used as input to plot the network using cytoscape
  needs to be run with:
    areaThreshold2.mat
      holds the interaction matrix of auxotroph pairs and their function annotation to group the strains
      
4. communityAnalysisM9.m 
  takes the RPM of members in assemblies grown over time in nutrient poor media (M9) and tracks and plots the composition over time. It also calculates the entropy (shannon diversity) of the assemblies over time
  needs to be run with:
    dataset_QS202.mat 
      holds the RPM of auxotrophs in 6 different assemblies
    functionAnnotations.xlsx 
      holds the names of auxotrophs used to start the assemblies and their function annotation
    colormap.mat
      holds a colormap for the different function annotations
      
5. communityAnalysisM9.m 
  takes the RPM of members in assemblies grown over time in nutrient rich media (LB) and tracks and plots the composition over time. It also calculates the entropy (shannon diversity) of the assemblies over time. Then, it plots them along the entropy of assemblies grown in nutrient poor media (M9)
  needs to be run with:
    RPM_COUNTS_QS10.mat
    functionAnnotations.xlsx
      same as in 4
    colormap.mat
      same as in 4
    entroM9.mat
      holds the entropy calculated in 4 for assemblies grown in nutrient poor media
      
6. plotCorrelationday14.m 
this script takes the correlation of assemblies' composition (grown in nutrient poor media (M9) or nutrient rich media (LB)) at day 14 and plots a heatmap of the correlation coefficients and a bar plot of the sorted correlation coefficients
  needs to be run with:
    corrLB.mat 
      holds the pearson correlation coefficient and the spearman correlation coefficient of auxotroph assemblies grown in nutrient rich media at day 14
    corrM9.mat
      holds the pearson correlation coefficient and the spearman correlation coefficient of auxotroph assemblies grown in nutrient poor media at day 14
      
7. randomSimulationSyntheticConsortia.m 
this script calculates an empirical p-value for getting an assembly supported by our pairwise interaction matrix (obtained in 2) given the 38 starting auxotroph strains
  needs to be run with:
    areaThreshold2.mat
      same as in 3
