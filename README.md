# ListeriaEx_scripts
Scripts used in our paper about genetic exchanges in Listeria's MGE.  
Only scripts used on the local dataset are provided, as the method is the same one on the global dataset.

## Citation
Please cite  "Genetic exchange networks bridge mobile DNA vehicles in the bacterial pathogen Listeria monocytogenes" by Muller et al. (Publication to come).

## MGEs annotation
 
Transposons were annotated with Bacant, in each assembly independently:
```
bacant -n assembly1.fa -o outputs_bacant_v3/assembly1/ -t 5 -p bacant-db-v3.1/ -d ResDB,IntegronDB,TransposonDB,repliconDB

#copy the output that sumarize everything, and add a colomn with the name of the assembly
awk -v var="assembly1" '{print $0, var}' OFS='\t' outputs_bacant_v3.1/assembly1/annotation.tsv > outputs_bacant_v3.1/assembly1/annotation_named.tsv

#Actually the output does not summarier everything: replicon in their on file
awk -v var="assembly1" '{print $0, var}' OFS='\t' outputs_bacant_v3.1/assembly1/replicon.tsv > outputs_bacant_v3.1/assembly1/replicon_named.tsv

#Once ran for all genomes, annotations were concatenated
cat  outputs_bacant_v3.1/assembly1/annotation_named.tsv  outputs_bacant_v3.1/assemblyN/annotation_named.tsv >  outputs_bacant_v3.1/all_annotation_named.tsv
```
Then, we used the script `extractTransposons.sh` to get the table of transposons only (`transposons.tsv`), and the fasta of transposons (`all_transposons.fasta`)


Plasmids were annotated with the MOB-recon tool from the MOB-suite, in each assembly independently.

Phage were annotated in all assemblies at once with the pipeline phageAnnotation, available at https://github.com/HeloiseMuller/phageAnnotation  
Two fasta of phage sequences were saved:  
* `phages_filtered_3TuningRemoval.fa` was filtered according to filters (i) and (ii) of the paper but not depending on the length.
* `phages_filtered_5000bp_AND_3TuningRemoval.fa` are the same phage sequences as above, but only keeping those which are at least 5000bp long.

Finally, all annotations were put together in a same file:  
`cat transposons/all_transposons.fasta plasmids/all_plasmids.fasta phages/phages_filtered_3TuningRemoval.fa > all_MGE.fasta`  
We also calculated the length of each MGE and saved the file as `all_MGE.length`.  

## Genes annotation
Genes were annotated with Prokka and AMRfinderplus

## Scripts to process the data
`ListeriaEx_network.R` finilizes the MGE annotation: it identifies PP elements and filter them out from the plasmid and phage annotations. It saves the final MGE annotation: `all_MGE_filtered.tsv`.  
It also generates the network in Figure 5 and Figure S9 of the paper.  
`ListeriaEx_connectivity.R` gives statistics on the network; it generates Figures S3 and S6 of the paper.  
`ListeriaEx_geneLocation.R` figures out the location of each gene: are their in any MGE?  
It also generates Figure 2 of the paper.  
`ListeriaEx_matrix.R` generates the matrix of Figure 1 of the paper. It also generates Figure 4 with the same data.  
`ListeriaEx_defenseSystem.R` generates Figure S1 of the paper. 
