Read Me


This is Kodee's MSc research project

I first started by filtering and subsetting my initial dataset using `workflow to cluster and subset.R`.
To help with clustering, I used `testing parameters to cluster COLD and CON datasets v2, Vim+ only.R`, `testing parameters on cluster1 + Ascl1+.R`, and then `testing gene markers on clean_combined_object.R`.

The sequencing metrics of the final dataset was performed using `sequencing metrics (for supplemental 2).R`.


The follow analysis / visualization was performed:
`cell cycle workflow.R`,
`expression umaps for thesis.R`,
`rna velocity workflow (done by SF).r`,
`pseudotime workflow.R`,
`augur workflow.R`, 
`cell proportion workflow.R`,
DEGs were determined for volcano plots and GO ORA using `deg for ORA workflow.R` and `deg analysis for thesis.R` was used to view shared DEG between treatment groups,
`volcano plot workflow v2.R` (just for Cycling NSPCs), `volcano plot workflow v2 (for other clusters).R` (for remaining clusters),
`go ora analysis.R`,
ligand-receptor analysis was done using `ligand receptor analysis (microglia data) workflow.R` to process the microglia data, then `cellchat individual clusters workflow.R` to run the analysis, and then `cellchat L-R Visualization_Formatted 2.R` for visualization. 


The control female verification dataset was processed and analyzed using:
`con f verification 2.R` to subset and cluster,
`testing gene markers on data with CON F verif dataset copy.R`,
`deg analysis for CON F Verification 2.R`,
`volcano plot for CON F Verification 2.R`.


GABAergic cells were processed and analyzed using:
`analyses of GABAergic cells.R`,
`testing markers on GABAergic cells.R`,
`expression maps for GABAergic cells.R`,
`gabaergic DEGs for ORA.R`,
`gabaergic DEGs comparisons.R`,
`gabaergic GO ORA Analysis.R`.


Glutamatergic cells were processed and analyzed using:
`analyses of glutamatergic cells.R`,
`testing markers on glutamatergic cells.R`,
`expression maps for glutamatergic cells.R`,
`glutamatergic DEGs for ORA.R`,
`glutamatergic DEGs comparisons.R`,
`glutamatergic GO ORA Analysis.R`.


Statistical analysis for neurosphere size distribution and RNAscope frequency curves was done using `art anova workflow 3.R`.

