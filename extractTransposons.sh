#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -v

#need seqtk activated

grep transposon  outputs_bacant_v3.1/all_annotation_named.tsv >  transposons.tsv 

while read i
do
        
	assembly=`echo $i | awk '{print $10}'`
        Tn=`echo $i | awk '{print $9}' | sed 's/.*gene=//'`

	#Make uniq ID for each transposon (I already checked that a same assembly never has twoce a Tn with the same name)
	ID=`paste -d : <(echo $assembly) <(echo $Tn)`

        #Write bed file
        echo $i |  awk -v var=$ID '{print $1,$4,$5,var}' OFS='\t'  > ${assembly}_${Tn}_transposon.bed

	#Extract fasta
        seqtk subseq  ${assembly}.fa ${assembly}_${Tn}_transposon.bed > ${assembly}_${Tn}_transposon.fasta

	#Does no care about my 4th column, so I have to add the ID myself
	sed -i "s/>.*/>$ID/" ${assembly}_${Tn}_transposon.fasta

	#Remove beds
	rm ${assembly}_${Tn}_transposon.bed

done <transposons.tsv

#Concatenate all results in one file
cat SRR*_transposon.fasta > all_transposons.fasta

#To be consistent with the other MGE fasta, I add ":transposon" at the end ID
sed -i '/^>/s/$/:transposon/'  all_transposons.fasta    
