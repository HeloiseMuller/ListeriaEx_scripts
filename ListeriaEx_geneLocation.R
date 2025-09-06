library(ape) #there is function to read gff
library(dplyr)
library(stringr)
library(data.table)
library(readxl)
library(ggplot)
library(ggpubr)
library(cowplot)

setwd("project_Listeria/")

####################################
##### Read all the annotations #####
####################################

##### 1st annotation = prokka

#List the 936 gff files without fasta
gff_files <- list.files("outputs_prokka", pattern = "_noseq.gff")

gffs <- lapply(gff_files, function(x){
    #read the gff
    dt <- read.gff(paste0("outputs_prokka/", x))
    #add some columns
    dt <- mutate(dt, 
        # add assembly name
        assembly = sub("_noseq.gff", "", x),
        # add gene name if there is one, otherwise NA
        gene = ifelse(str_detect(attributes, "gene="), sub(";.*", "", sub(".*gene=", "",attributes)), NA),
        # add production name
        product = ifelse(str_detect(attributes, "product="), sub(";.*", "", sub(".*product=", "", attributes)), NA),
        # also need ID, since it is in line with "Protein identifier" of AMRfinderplus
        ID = ifelse(str_detect(attributes, "ID="), sub(";.*", "", sub(".*ID=", "", attributes)), NA) 
        )
    #now that I exrtacted all the info I care about in $attribues, I remove this column (for conveniance purpose)
    #I also remove othet columns I don't use
    dt <- select(dt, -c(attributes, score, phase, source))

    return(dt)
    })

gffs <- rbindlist(gffs)

  
##### 2nd annotation = AMRfinderplus
# Use gff from prokka, and annotate genes of interest
# such as AMR, virulence, transposase (but not transposon) etc
amrfinder_files <- list.files("outputs_amrfinderplus", pattern = "_amrfinderplus.tsv")

amrfinder <- lapply(amrfinder_files, function(x){
    dt <- fread(paste0("outputs_amrfinderplus/", x))
    #space in colum name is annoying so I replace by _
    colnames(dt) <- gsub(" ", "_", colnames(dt))
    dt$assembly = sub("_amrfinderplus.tsv", "", x)
    return(dt)
})

amrfinder <- rbindlist(amrfinder)

##### 3rd annotation = bacant
# Also annotate AMR
# Plus transposons, integrons and replicon

bacant <- fread("all_annotation_named_DBv3.1.tsv", 
    col.names=c("contig", "element", "element2", "start", "end", "V6", "strand", "V8", "tag", "assembly"))  

#Note: it would be quite complicated to deal with AMR annotaton of that many tool,
#Bacant will thus be usefull to know where transposon and integron are.
#To note, bacant did not find a single integron in my dataset

##### 4th annotation = MOB-suites (done in script github_network936.R)
plasmid_filtered <- fread("plasmids_filtered_noPP.tsv")

##### 5th annotation = phageAnnotation (done in  github_network936.R)
phages_filtered <- fread("phages_filtered_noPP.tsv")

##########################
###### Deal with PP ###### 
########################## 

#Read MGE annotation, done in script github_network936.R)
allMGE <- fread("all_936SRR_MGE_filtered.tsv")
PP <- fread("PP_filtered.tsv")

#Get whether each contig is a plasmid or not :
plasmids_filtered <- allMGE[type=="plasmid" | type=="PP"] %>% mutate(., contig=sub(":plasmid", "", sub(".*_c", "c", MGEname)))
plasmids_filtered[, loc:=paste(assembly, contig, sep = ":")]


plasmids <- fread("correspondance_plasmids_mobilization.tsv",
                  header =F, col.names =c("sample_id", "size", "predicted_mobility"))

#split sample_id in assembly name & id plasmid
plasmids <- mutate(plasmids, assembly = sub(":.*", "", sample_id),
                   primary_cluster_id = sub(".*:", "", sample_id))


##########################
###### Get metadata ###### 
########################## 

##### Info about lineages and other group for each genome
meta <- data.table(read_excel("SupplementaryTables.xlsx", "Table_S1", skip=1))

###############################################

####################################
##### Put annotations together #####
####################################

#### 1st, resolve prokka & AMRfinderplus
# They have the same protein identifiers

#add suffix to all columns to see what comes from where
colnames(gffs) <- paste0(colnames(gffs)[1:length(gffs)], ".prokka")
colnames(amrfinder) <- paste0(colnames(amrfinder)[1:length(amrfinder)], ".amrfinder")

#join same gene of a same coordiantes
annot <- full_join(gffs, amrfinder, 
    by = c(
        "ID.prokka" = "Protein_identifier.amrfinder",
        "assembly.prokka" = "assembly.amrfinder",
        "seqid.prokka" = "Contig_id.amrfinder",
        "start.prokka" = "Start.amrfinder",
        "end.prokka" = "Stop.amrfinder",
        "strand.prokka" = "Strand.amrfinder"
        )
    ) %>% data.table


#rename columns which are in common between prokka & amrfinder
annot <- dplyr::rename(annot, contig = "seqid.prokka", start = "start.prokka", end = "end.prokka", assembly = "assembly.prokka", ID = "ID.prokka", strand = "strand.prokka")

#DOn't forget that 140 genes found by AMRfinderplu were not found by prokka
#So NA in all prokka columns
filter(annot, is.na(type.prokka))

#ANd much more genes were not found by AMRfinderplus since it only focuses on some

annot[, kind := 
        #if annotated by amrfinder, use its scope (just add "CDS ..."" before)
        ifelse(Scope.amrfinder!="", paste("CDS", Element_type.amrfinder), 
               ifelse(type.prokka=="CDS", "CDS others", type.prokka))]

#Here I favored amrfinder annotation
#But what of prokka genes with an AMR name not found by AMRfinder?

#Note: Element_type.amrfinder has 3 possible values: AMR, STRESS, or VIRULENCE
# Element_subtype.amrfinder gives a more detailed functional category
# Class.amrfinder :or AMR genes this is the class of drugs that this gene is known to contribute to resistance of.

annot[, loc:=paste(assembly, contig, sep = ":")]

#### 2nd Add info about molecule type (plasmid or chromosome)
annot$molecule_type = "chromosome"
annot[loc %in% plasmid_filtered$contig]$molecule_type = "plasmid"
annot[loc %in% PP$contig]$molecule_type = "PP"


#### 4th Are some genes in transposons?

#For more efficiency keep only columns and need
annot_l <- annot[, c("contig", "start", "end", "assembly")]

#and use integer because faster
#NOTE sub() is quite slow. 
#So it is much more efficient to do it only on unique values,
#and after to the correspondence
corresContig <- data.table(contig = unique(annot_l$contig), contigID = as.integer(sub("contig", "", unique(annot_l$contig))))
corresAssembly <- data.table(assembly = unique(annot_l$assembly), assemblyID = as.integer(sub("SRR", "", unique(annot_l$assembly))))
annot_l[, contig := corresContig[chmatch(annot_l$contig, contig), contigID]]
annot_l[, assembly := corresAssembly[chmatch(annot_l$assembly, assembly), assemblyID]]
#start & end are already integer

#Actually I could even do a unique ID for each comb of assembly-contig
annot_l[, ID := assembly*1000+contig]
#max value of contig is 153, so I need to add 3 0s at the end of assembly
#eg if assembly = 2782829 and contig : 1, we want ID = 2782829001

#These 2 columns are not useful anymore
annot_l <- annot_l[, -c("assembly", "contig")]

 
#Function that look whether overlap between 1 gene and a table of transposon
genes_inTn_f <- function(gene){
  
  #NOTE with apply
  #transposons[ID == gene[3],] will be empty no matter the values
  #so I firsly need to stock value of gene[3] (and. others) in variable
  
  gene_ID = gene["ID"]
  gene_start = gene["start"]
  gene_end = gene["end"]
  
  #candidate = same contig of the assembly as the gene we are looking at
  candidates <- transposons[ID == gene_ID,]
  
  #Check whether co localisation
  candidates[, loc := 
               #gene fully in the transposon:
               ifelse(start<=gene_start & end>=gene_end, "Full",  
                      #gene partially in transposon
                      ifelse((start>gene_start & start <gene_end) | #gene start before Transposon but finish inside
                               (end>gene_start & end<gene_end), "Overlap",  #gene start in Transposon but finish outside
                             #gene not in Transposon at all
                             "no")) 
  ]
  
  #Still candidate are Overlap or Full
  candidates <- candidates[loc!="no",]
  
  if(nrow(candidates)==0){
    out <- "FALSE"
  } else if(nrow(candidates)==1){
    out <- paste("TRUE", candidates$gene, candidates$loc, sep=":")
  } else {
    out <- paste("TRUE", "Tn_uncertain", ifelse(length(unique(candidates$log)==1), candidates$log, "NA"), sep=":")
  }
  
  cat("-")
  
  return(out)
}

#DO the same things as done with annot:
transposons <- bacant[element=="transposon",]
transposons[, gene:= splitToColumns(splitToColumns(transposons$tag, ";", column=2), "=", column=2)]

#Replace TnUncertain
transposons[assembly == "SRR1520067" & gene == "Tn925"]$gene = "TnUncertain"
transposons <- transposons[!(assembly == "SRR1520067" & gene!="TnUncertain")]

transposons[, contig := corresContig[chmatch(transposons$contig, contig), contigID]]
transposons[, assembly := corresAssembly[chmatch(transposons$assembly, assembly), assemblyID]]
transposons[, ID := assembly*1000+contig]
transposons <- transposons[, c("ID", "start", "end", "gene")]


res <- apply(X = annot_l, 1, genes_inTn_f)
res <- cbind(annot, inTransposon = res)

#hamonized gene names
#FOr now when 2 differente names, choose the one of amrfinder
res[, gene := ifelse(is.na(Gene_symbol.amrfinder), gene.prokka, Gene_symbol.amrfinder)]


--------------
  
#phages

#CHANGE THOS FUNCTON SO WORKS both wth transposon and pahes
  #NOT FINISHED
genes_inElement_f <- function(gene, element, table){
    #element can be "transposon" or "phage"
    #table is either the table tansposons or the table phages
    
    #NOTE with apply
    #transposons[ID == gene[3],] will be empty no matter the values
    #so I firsly need to stock value of gene[3] (and. others) in variable
    
    gene_ID = gene["ID"]
    gene_start = gene["start"]
    gene_end = gene["end"]
    
    #candidate = same contig of the assembly as the gene we are looking at
    candidates <- table[ID == gene_ID,]
    
    #Check whether co localisation
    candidates[, loc := 
                 #gene fully in the element:
                 ifelse(start<=gene_start & end>=gene_end, "Full",  
                        #gene partially in element
                        ifelse((start>gene_start & start <gene_end) | #gene start before element but finish inside
                                 (end>gene_start & end<gene_end), "Overlap",  #gene start in element but finish outside
                               #gene not in element at all
                               "no")) 
    ]
    
    #Still candidate are Overlap or Full
    candidates <- candidates[loc!="no",]
    
    if(nrow(candidates)==0){
      out <- "FALSE"
    } else if(nrow(candidates)==1){
      out <- paste("TRUE", candidates$gene, candidates$loc, sep=":")
    } else {
      out <- paste("TRUE", "uncertain", ifelse(length(unique(candidates$log)==1), candidates$log, "NA"), sep=":")
    }
    
    cat("-")
    
    return(out)
  }

#DO the same thing with phages
#WHAT start end should I use (many coordiantes in the table) ?
#https://github.com/jiarong/VirSorter2/issues/70 :
#"the trim_bp_start and trim_bp_end are the boundaries used for final viral contigs.
 #Be aware that host region trimming in VS2 is conservative, meaning there might be host region left.
  #You can use specialized prophage extraction tools to clean up."
phages_filtered_l <- phages_filtered[, c("assembly", "seqname", "contig", "start_inContig", "end_inContig")] %>% 
  #Rename this column so it is in adequation with transposon (we use the same function for both)
  dplyr::rename(., gene = "seqname", start = "start_inContig", end = "end_inContig") %>%
  mutate(., contig = sub(".*:", "", contig))
phages_filtered_l[, seqname := corresContig[chmatch(phages_filtered_l$contig, contig), contigID]]
phages_filtered_l[, assembly := corresAssembly[chmatch(phages_filtered_l$assembly, assembly), assemblyID]]
phages_filtered_l[, ID := assembly*1000+seqname]
phages_filtered_l <- phages_filtered_l[, -c("assembly", "seqname", "contig")]

res2 <- apply(X = annot_l, 1, genes_inElement_f, element = "phage", table = phages_filtered_l)
res2 <- cbind(res, inPhage = res2)

res2[inPhage!=F,]

#Not a single gene both in transposon & phage
res2[inTransposon!=F & inPhage!=F,]

#So we can have just one column for both:
res2 <- mutate(res2, element = ifelse(inTransposon==F & inPhage==F, "no", ifelse(inTransposon!=F, "transposon", "phage")))

#Since never both in transpon & phage, I can summarize this info in a single column
res2[inTransposon!=F & inPhage!=F,] #indeed, empty
res2[, MGE := ifelse(inPhage != F, "phage", ifelse(inTransposon!=F, "transposon", NA))]


###########
##Prokka never give twice the same gene name
#So if 2 copies of metE in a genome, one will be called metE_1, and the other met_2
#But if another genomes has only 1 vopy of metEn it will be called metE
#--> number metE_1 <= metE_2 <= metE_3, etc but metE can have any value

res2[, geneUndup:=sub("_.*", "", gene)]


#Add info about whether this gene is a gene duplicated in THIS assembly
res2[, nbCopy_gene.prokka := ifelse(str_detect(gene, "_"), "duplicated", ifelse(gene.prokka=="", NA, "single"))]

       
annot <- res2

################
# SOme gnes names and categories have to be corrected
##############
annot[geneUndup== "bcrC" | geneUndup=="bcrC"] #actually no need to modify this as I don't have figures with Bacitracin. Only modify manuscript

annot[geneUndup=="vga(G)"]$geneUndup = "abc-f"

annot[geneUndup=="emrC"]$kind = "CDS STRESS"
annot[geneUndup=="emrC"]$Element_type.amrfinder = "STRESS"
annot[geneUndup=="emrC"]$Element_subtype.amrfinder = "BIOCIDE"

fwrite(annot, "harmionizedAnnotation_correctGeneNames.tsv", sep='\t')

annotCare <-  annot[gene!="blaZ",] #reviewer does not want blaZ

###############################

###############################
########### FIGURE 2 ##########
####### of the article ####### 
###############################

### PANEL A

annot_summarize <- annotCare[Element_type.amrfinder!=""] %>%
  #Group by assembly and scope
  group_by(phylogenetic_lineage, Element_type.amrfinder, geneUndup) %>%
  #count the number of genes for each group
  dplyr::summarize(nb=length(unique(assembly))) %>% 
  data.table() 
  #When 0 gene for a group, no line. We need to add one

dummy <- data.table(Element_type.amrfinder = rep(distinct(annotCare[Element_type.amrfinder!="", c("Element_type.amrfinder", "geneUndup")])$Element_type.amrfinder , each = 4),
           geneUndup = rep(distinct(annotCare[Element_type.amrfinder!="", c("Element_type.amrfinder", "geneUndup")])$geneUndup , each = 4),
           phylogenetic_lineage = rep(unique(meta$phylogenetic_lineage), nrow(distinct(annotCare[Element_type.amrfinder!="", c("Element_type.amrfinder", "geneUndup")])))
)

annot_summarize <- left_join(dummy, annot_summarize) %>%
  mutate(., nb = ifelse(is.na(nb), 0, nb)) %>%
  #Calculate percentage of lieage with each gene
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage"="V1")) %>%
  mutate(., perc_ofLineage = nb*100/N)

f2a <- ggplot(annot_summarize, aes(x = geneUndup, y = perc_ofLineage, fill = phylogenetic_lineage)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~Element_type.amrfinder, scale = "free", space = "free") + 
  theme_bw() + 
  theme(axis.text=element_text(angle=90, size=8)) + 
  ylim(0,100) +
  ylab("Percentage of assemlbly") +
  xlab("Genes") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442", "#CC79A7"), name = "Lineage")

### PANEL B
f2b <- allMGE[type=="phage"] %>% group_by(phylogenetic_lineage) %>%
  summarize(nAssembly = length(unique(assembly))) %>% 
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage" = "V1")) %>%
  mutate(., perc = nAssembly*100/N) %>% ggplot(., aes(x=phylogenetic_lineage, y = perc, fill = phylogenetic_lineage))  +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  ylab("Percentage of assembly") +
  xlab("Phylogenetic lineages") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442", "#CC79A7"), name = "Lineage") +
  theme(legend.position = "none")

### PANEL C

f2c <- allMGE[type=="plasmid"] %>% 
  mutate(., name = splitToColumns(nodeName, split = ":", columns = 2)) %>%
  #Count the number of assembly per lineaage with each plasmid
  group_by(name, phylogenetic_lineage) %>%
  summarize(nAssembly = length(unique(assembly))) %>% 
  
  #Add a column saying how many assemblies in total for that lineage
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage" = "V1")) %>%
  #Calculate percentage of assemblies in this lineage with that plasmid
  mutate(., perc = nAssembly*100/N) %>%
  
  #Make the figure
  ggplot(., aes(x=reorder(name, -nAssembly), y = perc, fill = phylogenetic_lineage)) + #y can be number of nAssembly or perc
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + 
  ylab("Percentage of assembly") +
  xlab("Plasmid group (n = 21)")  +
  scale_fill_manual(values = c("#56B4E9", "#009E73"), name = "Lineage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

### PANEL D

f2d <- allMGE[type=="transposon"] %>% 
  mutate(., name = splitToColumns(nodeName, split = ":", columns = 2)) %>%
  group_by(name, phylogenetic_lineage) %>% 
  summarize(nAssembly = length(unique(assembly))) %>% 
  
  #Add a column saying how many assemblies in total for that lineage
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage" = "V1")) %>%
  #Calculate percentage of assemblies in this lineage with that plasmid
  mutate(., perc = nAssembly*100/N) %>%
  
  ggplot(., aes(x=name, y = perc, fill = phylogenetic_lineage)) + #y can by perc or nAssembly
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + 
  ylab("Percentage of assembly") + xlab("Transposon family (n = 6)")  +
  scale_fill_manual(values = c("#56B4E9", "#009E73"), name = "Lineage") +
  theme(legend.position = "none")

#### Put figure together
gt <- arrangeGrob(f2a, f2b, f2c, f2d, 
                  ncol = 2, nrow=2, 
                  layout_matrix = rbind(c(1,1,1,1,2), c(3,3,3,4,4)))

pdf("Figure2.pdf", width = 12, height = 10)
as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.8, 0, 0.6), y = c(1, 1, 0.5, 0.5)) # Add labels
dev.off()
