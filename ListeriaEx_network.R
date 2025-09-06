library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(igraph)
library(circlize)


setwd("project_Listeria/")

meta <- data.table(read_excel("SupplementaryTables.xlsx", "Table_S1", skip=1))

#################
### MGE annotation

#Read phages annotation
#THis is the 1st output of the second script of phageAnnotation.R
#Those phages were not filtered but they contains all information needed for filtering
phages <- fread("all_936SRR_phages_quality.tsv", sep = "\t")

#Read the file containing length of all putative MGE
#This file also contains phages that we filtered out (tot phage before filtering = 2728) and the plasmid AG939
#This file was done before PP annotation
allMGE <- fread("all_936SRR_MGE.length", header = F, col.names = c("MGEname", "length"))
allMGE[, type:=splitToColumns(MGEname, split = ":", columns = 3)]
#Tot MGE here is 3554

#According to side analyses, the huge "plasmid" AG939 is probably a false positive
allMGE <- allMGE[!str_detect(MGEname, "AG939")] #We now have 437 contigs annotated as plasmids

#Extract plasmid annotation from allGME
plasmid_filtered <- filter(allMGE, type=="plasmid") %>% mutate(., contig=sub(":.*_", ":", str_extract(MGEname, "([^:]*:[^:]*)")))
plasmid_filtered[, assembly := sub(":.*", "", MGEname)]
plasmid_filtered[, tmp := splitToColumns(MGEname, split = ":", columns = 2)]
plasmid_filtered[, grp := sub("_contig.*", "", tmp)]
plasmid_filtered[, reconstruct := paste(assembly, grp, sep = ":")]
plasmid_filtered <- plasmid_filtered[, -"tmp"]
plasmid_filtered <- left_join(plasmid_filtered, meta[,c("ID", "phylogenetic_lineage")], by = c("assembly"="ID"))

                        
#Replace contigs that could be phage or Pl by PP
allMGE[MGEname %in% plasmid_filtered[contig %in% phages[length==length.plasmid]$contig]$MGEname]$type = "PP" #Used plasmid annotation for PP annotation
allMGE <- allMGE[!MGEname %in% paste0(phages[length==length.plasmid]$contigVir, ':phage')] #Remove phage annotated as PP from annotation
#I could have done the opposite: using phage annotation for PP annotation and then remove plasmid annotation 

#Now that PP are annotated, only keep phage passing our filters (2nd output of the pipeline phageAnnotation)
phages_filtered <- fread("phages_filtered_5000bp_AND_3TuningRemoval.tsv")
allMGE <- allMGE[type!="phage" | MGEname %in% paste0(phages_filtered$seqname , ":phage")]


#Rename TnUncertain
allMGE[MGEname=="SRR1520067:Tn925:transposon",]$MGEname = "SRR1520067:TnUncertain:transposon"
allMGE <- allMGE[!str_detect(MGEname, "SRR1520067") | (str_detect(MGEname, "phage") | str_detect(MGEname, "TnUncertain"))]


#Rename "novel" plasmids
allMGE[MGEname=="SRR14524738:novel_6d6c6c02c6d2451d503888c6577b537c_contig00026:plasmid"]$nodeName = "SRR14524738:novel_A:plasmid"
allMGE[MGEname=="SRR14404492:novel_c2ff2894c49bc57808750dee62b74edf_contig00021:plasmid"]$nodeName = "SRR14404492:novel_B:plasmid"
allMGE[MGEname=="SRR14404492:novel_c2ff2894c49bc57808750dee62b74edf_contig00022:plasmid"]$nodeName = "SRR14404492:novel_B:plasmid"
allMGE[MGEname=="SRR13744454:novel_ce9e32fd426b728011b1a6ef96f50b16_contig00017:plasmid"]$nodeName = "SRR13744454:novel_C:plasmid"

##############################
#### Add info in all MGE #####
## And save MGE annotations ##
##############################

allMGE$nodeName = allMGE$MGEname
allMGE[type=="plasmid", ]$nodeName = sub("_contig.*", ":plasmid" , allMGE[type=="plasmid", ]$nodeName)

length(unique(allMGE$nodeName)) #2332 nodes

#Add assembly and lineage
allMGE[, assembly:=sub(":.*", "", MGEname)]
allMGE <- left_join(allMGE, meta[, c("ID", "phylogenetic_lineage")], by = c("assembly"="ID"))

#Set colors
set_colLinege <- data.table(lineage = c("I", "II", "III", "IV"), col = c("#56B4E9", "#009E73","#F0E442",  "#CC79A7"))
allMGE[, colLinege:=set_colLinege[match(phylogenetic_lineage, lineage), col,]]
allMGE[, colType:=ifelse(type=="transposon", "blue", ifelse(type=="phage", "red", ifelse(type=="PP", "orange", "gray")))]

fwrite(allMGE, "all_936SRR_MGE_filtered.tsv", sep = "\t")

#Also save filtered phage and plasmids for the annotation
PP <-  plasmid_filtered[MGEname %in% allMGE[type=="PP"]$MGEname]
plasmid_filtered <- plasmid_filtered[MGEname %in% allMGE[type=="plasmid"]$MGEname]
fwrite(plasmid_filtered, "plasmids_filtered_noPP.tsv", sep = "\t")
fwrite(PP, "PP_filtered.tsv", sep = "\t")
fwrite(phages_filtered, "phages_filtered_noPP.tsv")


##############################
#### Blast that compare #####
##### each pair of MGE ######
##############################

#THe blast was done before filtered MGE, in case decide to change filters later
#THis is why in this blast we have 3554 unique c(query, subject)
blastMGE <- fread("all_936SRR_MGE_selfBlastn.out", 
                  col.names = c("query", "subject", "pident", "length",  "qlen", "slen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))



blastMGE <- blastMGE[query!=subject,]
blastMGE_long <- blastMGE[length>=200,]

#Put reciprocal in the same order
blastMGE_long <- mutate(blastMGE_long, 
                        Q = ifelse(query<subject, query, subject),
                        S =  ifelse(query<subject, subject, query),
                        Qlen =  ifelse(query<subject, qlen, slen),
                        Slen =  ifelse(query<subject, slen, qlen),
                        Qstart =  ifelse(query<subject, qstart, sstart),
                        Qend =  ifelse(query<subject, qend, send),
                        Sstart =  ifelse(query<subject, sstart, qstart),
                        Send =  ifelse(query<subject, send, qend)
) %>%
  select(., -c("query", "subject", "qlen", "slen", "qstart", "qend", "sstart", "send")) %>%
  rename(., query = "Q", subject = "S")



#Get assembly name & MGE name & MGE type
blastMGE_long[, Qassembly :=  splitToColumns(query, ":", columns = 1)]
blastMGE_long[, Sassembly :=  splitToColumns(subject, ":", columns = 1)]

# There is a problem with transposons of SRR1520067: say 6 put it is actually the same coordinates --> TnUncertain
blastMGE_long[Qassembly == "SRR1520067" & str_detect(query, "transposon"),]$query = "SRR1520067:TnUncertain:transposon"
blastMGE_long[Sassembly == "SRR1520067" & str_detect(query, "transposon"),]$subject = "SRR1520067:TnUncertain:transposon"

#At this stage, we only keep MGE we retained in all_936SRR_MGE_filtered.tsv
blastMGE_long <- blastMGE_long[query %chin% allMGE$MGEname & subject %chin% allMGE$MGEname]
#It deletes AG939, phage of bad quality and redundance due to PP

#Add info about type each MGE
blastMGE_long  <- left_join(blastMGE_long, allMGE[, c("MGEname", "type")], by = c("query"="MGEname")) %>%
  left_join(.,  allMGE[, c("MGEname", "type")], by = c("subject"="MGEname"))  %>%
  rename(., Qtype="type.x", Stype="type.y")

#Add info about kind of interaction (eg phage-phage or plasmid-phage)
blastMGE_long[, pairType:= ifelse(Qtype<Stype, paste0(Qtype, "-", Stype), paste0(Stype, "-", Qtype))]

#Give same names to same MGE family 
#For phage, since were not clusterized I need an individual name, so I have to keep assembly name
#For plasmid, get rid of contig name so we can keep just one per plasmid cluster
blastMGE_long[, QMGEname := ifelse(Qtype=="PP", sub(":.*_", ":", sub(":plasmid", "", query)),
                                   ifelse(Qtype=="phage",  sub(":phage", "", query),
                                          ifelse(Qtype=="plasmid", sub("_.*", "", splitToColumns(query, ":", columns = 2)),
                                                 ifelse(Qtype=="transposon", splitToColumns(query, ":", columns = 2), NA))))]
blastMGE_long[, SMGEname := ifelse(Stype=="PP", sub(":.*_", ":", sub(":plasmid", "", subject)),
                                   ifelse(Stype=="phage",  sub(":phage", "", subject),
                                          ifelse(Stype=="plasmid", sub("_.*", "", splitToColumns(subject, ":", columns = 2)),
                                                 ifelse(Stype=="transposon", splitToColumns(subject, ":", columns = 2), NA))))]



#When select best hit per pair, we want to take into account that a same plasmid is split on several contigs
blastMGE_long <- mutate(blastMGE_long, queryInitial = query, subjectInitial = subject)
blastMGE_long[Qtype=="plasmid", ]$query = sub("_contig.*", ":plasmid" , blastMGE_long[Qtype=="plasmid", ]$query)
blastMGE_long[Stype=="plasmid", ]$subject = sub("_contig.*", ":plasmid" , blastMGE_long[Stype=="plasmid", ]$subject)

#All this handling led to some query=subject. Delete them
blastMGE_long <- blastMGE_long[query != subject,]

#order by best score
setorder(blastMGE_long, -bitscore)

#Keep the best score per query-subject
blastMGE_best <- blastMGE_long[!duplicated(paste(query, subject))]


#Rename "novel" plasmid
blastMGE_best$query <- sub("6d6c6c02c6d2451d503888c6577b537c", "A", blastMGE_best$query)
blastMGE_best$query <- sub("c2ff2894c49bc57808750dee62b74edf", "B", blastMGE_best$query)
blastMGE_best$query <- sub("ce9e32fd426b728011b1a6ef96f50b16", "C", blastMGE_best$query)

blastMGE_best$subject <- sub("6d6c6c02c6d2451d503888c6577b537c", "A", blastMGE_best$subject)
blastMGE_best$subject <- sub("c2ff2894c49bc57808750dee62b74edf", "B", blastMGE_best$subject)
blastMGE_best$subject <- sub("ce9e32fd426b728011b1a6ef96f50b16", "C", blastMGE_best$subject)

#Save this for further use
fwrite(blastMGE_best, "all_936SRR_MGE_selfBlastn_best.tbl", sep = "\t")

##############################
##### Build the network  ##### 
##############################

#allMGE has one line per contigs. We want 1 node per reconstructed plasmid
nodes <- allMGE[!duplicated(allMGE$nodeName)]
#Put column in right order and keep only columns of interest (don't kee length, as does not maje sense anymore for plasmid!)
nodes <- nodes[, c("nodeName", "type", "assembly", "phylogenetic_lineage", "colLinege", "colType")]


net <- graph_from_data_frame(d = blastMGE_best[, c("query", "subject", "bitscore", "pident", "pairType")], 
                             vertices = nodes, 
                             directed=F)

#Get clusters
cfg <- cluster_fast_greedy(net, weights = log10(E(net)$bitscore))

#Look who is in what community
cfg_summary <- data.frame(membership = cfg$membership, names = cfg$names) %>% group_by(membership, names) %>%summarise(n=n()) %>% data.table()

#Is there any MGE with the same name not in a same com?
cfg_summary %>% 
  group_by(names) %>%
  summarize(nMembership = length(unique(membership))) %>% 
  filter(., nMembership>1) #a same name is alwys in the same com

#I can also add membership in table allMGE
nodes[, membership:=cfg_summary[match(nodeName, names), membership,]]

fwrite(nodes, "all_936SRR_MGE_filtered_nodes.tbl", sep="\t")

#Run the network again with info about memberships
net <- graph_from_data_frame(d = blastMGE_best[, c("query", "subject", "bitscore", "pident", "pairType")], 
                             vertices = nodes, 
                             directed=F)


#Set width edges
#max bitscore is 15697, which is way too bit for plot
#log (instead of division) will decrease diff between min and max
E(net)$weight <-  log10(E(net)$bitscore)


#Save my graph so I can read it on  Cytoscape
write_graph(graph = net, file = "graph_all_936SRR_MGE_filtered.gml", format = "gml")


#########################################################
#########################################################

###########################
## Any link between CC and network?
# This part will generate the figure S9
###########################

#Add info about CC
blastMGE_best <- left_join(blastMGE_best, meta[, c("ID", "clonal_complex")], by = c("Qassembly"="ID")) %>%
  left_join(., meta[, c("ID", "clonal_complex")], by = c("Sassembly"="ID"), suffix = c(".q", ".s")) %>%
  mutate(., pairCC = ifelse(clonal_complex.q<clonal_complex.s, paste(clonal_complex.q, clonal_complex.s, sep = "-"),  paste(clonal_complex.s, clonal_complex.q, sep = "-")))

allMGE <- left_join(allMGE, meta[, c("ID", "clonal_complex")], by = c("assembly"="ID"))

MGE_perCC <- allMGE %>% group_by(clonal_complex, type) %>% summarise(nMGE = n()) %>% data.table() %>%
  mutate(., CC_type = paste(clonal_complex, type, sep = "-"))

compMGE_perCC <- data.table(CC_type_q = rep(MGE_perCC$CC_type, each = nrow(MGE_perCC)), 
                            nMGE_q = rep(MGE_perCC$nMGE, each = nrow(MGE_perCC)))
compMGE_perCC$CC_type_s = ""
compMGE_perCC$nMGE_s = 0

for(i in 1:nrow(MGE_perCC)){
  for (j in 1:nrow(MGE_perCC)){
    #print(paste0("i= ",i))
   # print(paste0("j= ",j))
    
    compMGE_perCC[CC_type_q==MGE_perCC[i,]$CC_type][j,]$CC_type_s = MGE_perCC[j,]$CC_type
    compMGE_perCC[CC_type_q==MGE_perCC[i,]$CC_type][j,]$nMGE_s = MGE_perCC[j,]$nMGE
  }
}

compMGE_perCC <- mutate(compMGE_perCC, pair = ifelse(CC_type_q<CC_type_s, paste(CC_type_q, CC_type_s, sep = ":"), paste(CC_type_s, CC_type_q, sep = ":"))) %>%
  mutate(., comb = nMGE_q*nMGE_s) %>%
  filter(., CC_type_q!=CC_type_s)
compMGE_perCC <- compMGE_perCC[!duplicated(compMGE_perCC$pair)]

blastMGE_perCC <- blastMGE_best %>% filter(.,  Stype!=Qtype & pairType!="plasmid-transposon") %>% 
  mutate(., CC_type_q = paste(clonal_complex.q, Qtype, sep = "-"), CC_type_s = paste(clonal_complex.s, Stype, sep = "-")) %>%
  mutate(., pair = ifelse(CC_type_q<CC_type_s, paste(CC_type_q, CC_type_s, sep = ":"), paste(CC_type_s, CC_type_q, sep = ":"))) %>%
  mutate(., virulence_q = ifelse(clonal_complex.q %in% c("CC1",  "CC6", "CC5", "CC4", "CC217", "CC2"), "hyper", ifelse(clonal_complex.q == "CC121", "hypo", "other")), 
         virulence_s = ifelse(clonal_complex.s %in% c("CC1",  "CC6", "CC5", "CC4", "CC217", "CC2"), "hyper", ifelse(clonal_complex.s == "CC121", "hypo", "other"))) %>%
  mutate(., pairVirulence = ifelse(virulence_q<virulence_s, paste(virulence_q, virulence_s, sep = ":"), paste(virulence_s, virulence_q, sep = ":"))) %>%
  group_by(pairCC, pairType, pairVirulence, pair) %>% summarise(nCo=n()) %>% data.table() 

blastMGE_perCC <- left_join(blastMGE_perCC, compMGE_perCC[, c("pair", "comb")], by = "pair") %>%
  mutate(., perc_co = nCo*100/comb)
  

#Identify a CC with more co
blastMGE_perCC$focus1 = splitToColumns(blastMGE_perCC$pairCC, split = "-", columns = 1)
blastMGE_perCC$focus2 = splitToColumns(blastMGE_perCC$pairCC, split = "-", columns = 2)
blastMGE_perCC[, compCC := ifelse(focus1==focus2, "intraCC", "interCC")] 
blastMGE_perCC[focus1==focus2 & pairType=="phage-plasmid"] #Only 3 CC have such exchanges intra CC
blastMGE_perCC_focus <- blastMGE_perCC[, -c("nCo", "comb", "focus2")] %>% 
  rename(., focus2="focus1") %>%
  rbind(., blastMGE_perCC[, -c("nCo", "comb", "focus1")]) %>%
  rename(., focus="focus2")

#intra CC are uselessly duplicated
blastMGE_perCC_focus <- distinct(blastMGE_perCC_focus)

#Give colors to CC in accordance with rest of the paper
colCircos_rightformat <- c(CC1180 = "#CCCCCC", CC89 = "#CCCCCC",  ST1064 = "#636161", CC121 = "#636161", CC11 = "#D53E4F", CC1907 = "#CCCCCC", CC155 = "#636161",
                           ST376 = "#636161", CC207 = "#636161", CC7 = "#636161", CC321 = "#636161", CC364 = "#CCCCCC", CC8 = "#CCCCCC",CC671 = "#636161",
                           CC415 = "#CCCCCC",CC635 = "#CCCCCC", CC9 = "#636161", CC224 = "#CCCCCC", CC288 = "#636161", ST2038 = "#CCCCCC",
                           CC5 = "#66C2A5", CC3 = "#CCCCCC",   CC736 = "#636161",CC59 = "#CCCCCC", CC379 = "#636161",CC87 = "#CCCCCC", CC88 = "#636161",
                           CC4 = "#ABDDA4", CC2 = "#F46D43", CC217 = "#FDAE61", CC6 = "#5E4FA2",  CC1 = "#9E0142")


SF9_a <- ggplot(blastMGE_perCC_focus[pairType=="phage-plasmid"], aes(x=focus, y=perc_co)) +
  geom_boxplot() + geom_point(aes(col = pairVirulence), cex = 2)+
  ylab("Percentage of possible phage-plasmid connections") +
  theme_bw() +
  xlab("Focal CC") +
  #Make the big figure without legend for scale reason, I'll add it with Inkscape
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

SF9_b <- ggplot(blastMGE_perCC_focus[pairType=="phage-plasmid"], aes(x=focus, fill = focus)) + geom_bar(stat = "count") +  
  ylab("Number of CCs") +
  xlab("Focal CC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = colCircos_rightformat) +
  theme(legend.position = "none") 


## Make a chord diagram for panel C

#Generate some data to be able to build this chord diagral
allMGE_PhPl <- filter(allMGE, type=="plasmid" | type=="phage") %>% 
  mutate(., ID_assembly= paste(type, clonal_complex, assembly, sep = "-")) %>% 
  mutate(., ID= paste(type, clonal_complex, sep = "-")) %>% 
  select(., c("ID", "ID_assembly")) %>%
  distinct() %>%
  group_by(ID) %>% summarise(nAssembly = n()) %>%
  data.table()


pairsID <- data.table(from = rep(allMGE_PhPl[str_detect(ID, "phage")]$ID, each = nrow(allMGE_PhPl[str_detect(ID, "plasmid")])), 
                      from_nAssembly = rep(allMGE_PhPl[str_detect(ID, "phage")]$nAssembly, each = nrow(allMGE_PhPl[str_detect(ID, "plasmid")])),
                      to ="",
                      to_nAssembly=0)
for(i in 1:nrow(allMGE_PhPl[str_detect(ID, "phage")])){
  for(j in 1:nrow(allMGE_PhPl[str_detect(ID, "plasmid")])){
    pairsID[from==allMGE_PhPl[str_detect(ID, "phage")][i,]$ID][j,]$to = allMGE_PhPl[str_detect(ID, "plasmid")][j,]$ID
    pairsID[from==allMGE_PhPl[str_detect(ID, "phage")][i,]$ID][j,]$to_nAssembly = allMGE_PhPl[str_detect(ID, "plasmid")][j,]$nAssembly
    
  }
}
pairsID[, pairAssembly_max := from_nAssembly*to_nAssembly]
pairsID[, pair:=paste(from, to, sep = ":")]

pairsID[, fromCC:=sub(".*-", "", from)]
pairsID[, toCC:=sub(".*-", "", to)]
pairsID[, pairCC:=ifelse(fromCC<toCC, paste(fromCC, toCC, sep = "-"),  paste(toCC, fromCC, sep = "-"))]

pairsID <- pairsID %>% group_by(pairCC) %>% summarise(pairAssembly_max=sum(pairAssembly_max)) %>% data.table()

#Build the chord diagram
forCircos <- blastMGE_best[pairType=="phage-plasmid",] %>% 
  mutate(., CC_type_q = paste(Qtype, clonal_complex.q, sep = "-"), CC_type_s = paste(Stype, clonal_complex.s, sep = "-")) %>%
  mutate(., pair = ifelse(CC_type_q<CC_type_s, paste(CC_type_q, CC_type_s, sep = ":"), paste(CC_type_s, CC_type_q, sep = ":"))) %>%
  mutate(., virulence_q = ifelse(clonal_complex.q %in% c("CC1",  "CC6", "CC5", "CC4", "CC217", "CC2"), "hyper", ifelse(clonal_complex.q == "CC121", "hypo", "other")), 
         virulence_s = ifelse(clonal_complex.s %in% c("CC1",  "CC6", "CC5", "CC4", "CC217", "CC2"), "hyper", ifelse(clonal_complex.s == "CC121", "hypo", "other"))) %>%
  mutate(., pairVirulence = ifelse(virulence_q<virulence_s, paste(virulence_q, virulence_s, sep = ":"), paste(virulence_s, virulence_q, sep = ":"))) %>%
  mutate(., pairAssembly = ifelse(Qassembly<Sassembly, paste(Qassembly, Sassembly, sep = ":"),  paste(Sassembly, Qassembly, sep = ":"))) %>%
  group_by(pairCC, pairType, pairVirulence, pair) %>% summarise(n_pairsAssemblies=length(unique(pairAssembly))) %>%
  group_by(pairCC, pairType, pairVirulence) %>% summarise(n_pairsAssemblies_tot=sum(n_pairsAssemblies)) %>% data.table()
                                                                  
forCircos <- left_join(forCircos, pairsID) %>%
  mutate(., perc_pairAssemblies = n_pairsAssemblies_tot*100/pairAssembly_max) %>%
  mutate(., from  = sub("-.*", "", pairCC), to = sub(".*-", "", pairCC))


#### Save figure S9


#I would have like to generate all 3 panels at once, but it does not work with circos
gt <- arrangeGrob(SF9_a, SF9_b, 
                  layout_matrix = rbind(c(1,1,1), c(1,1,1),c(2,2,2)))


pdf("SuppFig9_AB.pdf", width = 10, height = 8)

as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.35))
dev.off()

#Then save panel C appart. And put A,B,C together in inskape
pdf("SuppFig9_C.pdf", width = 5, height = 5)
chordDiagram(forCircos[, c("from", "to", "perc_pairAssemblies")], annotationTrack = c("name", "grid"),
             grid.col = colCircos_rightformat,
             transparency = 0.3,
             link.lwd = 0.3,    # Line width
             link.lty = 1,    # Line type
             link.border = "black",
             grid.border = "black")
dev.off()

