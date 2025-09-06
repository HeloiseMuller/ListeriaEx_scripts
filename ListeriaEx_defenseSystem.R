library(data.table)
library(dplyr)
library(readxl)

setwd("project_Listeria/")

annot <- fread("harmionizedAnnotation_correctGeneNames.tsv")

##Metadata
meta <- data.table(read_excel("SupplementaryTables.xlsx", "Table_S1", skip=1))

annot <- left_join(annot, meta[, c("ID", "phylogenetic_lineage")], by =c("assembly"="ID"))

annot_restrict <- annot[str_detect(product.prokka, "restriction")]  %>% 
  group_by(product.prokka, geneUndup, phylogenetic_lineage) %>%
  summarize(nAssembly=length(unique(assembly))) %>% data.table() %>%
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage"="V1")) %>%
  mutate(., perc_ofLineage = nAssembly*100/N) %>%
  filter(., geneUndup %in% c("mrr", "hsdR", "sau3AIR")) %>%
  mutate(., system = "Restriction-modification")

annot_crispr <- annot[str_detect(product.prokka, "CRISPR")] %>% 
  group_by(product.prokka, geneUndup, phylogenetic_lineage) %>%
  summarize(nAssembly=length(unique(assembly))) %>% data.table() %>%
  left_join(., data.table(table(meta$phylogenetic_lineage)), by = c("phylogenetic_lineage"="V1")) %>%
  mutate(., perc_ofLineage = nAssembly*100/N)  %>%
  mutate(., system = "CRISPR")


#Figure 
d1 <- rbind(annot_restrict, annot_crispr)  %>%
  filter(., phylogenetic_lineage!="IV") %>%
  ggplot(., aes(x=geneUndup, y = perc_ofLineage, fill = phylogenetic_lineage)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(~system,  scale = "free", space = "free") +
  theme_bw() + 
  ylim(0,100) +
  ylab("Percentage of assemlbly") +
  xlab("Genes associated with defense systems") +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442"), name = "Lineage")


#Phage
phages <- fread("annotVirus_all936_clean/phages_filtered.tsv")

d2 <- left_join(phages, meta[, c("ID", "phylogenetic_lineage")], by = c("assembly"="ID")) %>%
  filter(.,  phylogenetic_lineage!="IV") %>%
  group_by(assembly, phylogenetic_lineage) %>%
  summarise(sumLength = sum(length)) %>%
  ggplot(., aes(y=sumLength, x=phylogenetic_lineage, fill=phylogenetic_lineage)) +
  geom_boxplot() + ylab("Total length of phage in each assembly (bp)") + 
  theme_bw() + 
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#F0E442"), name = "") +
  xlab("Phylogenetic lineage" ) +
  guides(fill="none")


#### Put figure together
gt <- arrangeGrob(d1, d2, 
                  ncol = 2, nrow=1, 
                  layout_matrix = rbind(c(1,1,1,2)))

pdf("FigureS1.pdf", width = 10, height = 6)
as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.75), y = c(1, 1)) # Add labels
dev.off()


#### Statistical tests

## For phages:

phages_stats <-  left_join(phages, meta[, c("ID", "phylogenetic_lineage")], by = c("assembly"="ID")) %>%
  filter(.,  phylogenetic_lineage!="IV") %>%
  group_by(assembly, phylogenetic_lineage) %>%
  summarise(sumLength = sum(length)) %>%
  mutate(sumLength_log = log(sumLength))

pairwise.wilcox.test(phages_stats$sumLength_log, phages_stats$phylogenetic_lineage,
                     p.adjust.method = "BH")

## For defense system

restric_stats <-  annot[str_detect(product.prokka, "restriction") | str_detect(product.prokka, "CRISPR")]  %>%
  filter(.,  phylogenetic_lineage!="IV") %>%
  mutate(., system = ifelse(str_detect(product.prokka, "restriction"), "Restriction-modification", "CRISPR")) %>%
  group_by(assembly, system, phylogenetic_lineage) %>%
  summarize(nGene = length(unique(geneUndup))) %>% data.table() 

#Add assemblies with value of 0
restric_dummy <- data.table(meta[!ID %in% restric_stats[system=="Restriction-modification"]$assembly & phylogenetic_lineage!="IV", c("ID", "phylogenetic_lineage")],
                            system="Restriction-modification",
                            nGene = 0) %>%
  rbind(.,  data.table(meta[!ID %in% restric_stats[system=="CRISPR"]$assembly & phylogenetic_lineage!="IV", c("ID", "phylogenetic_lineage")],
                             system="CRISPR",
                             nGene = 0)
  ) %>%
  rename(., assembly="ID")

restric_stats <- rbind(restric_stats, restric_dummy)

restric_stats[, presence:=ifelse(nGene==0, F, T)]

#All categories have effective >5
restric_stats_summary <- restric_stats %>% group_by(phylogenetic_lineage, presence, system) %>%summarise(nEffectif = n()) %>% data.table()

#Run a chi2 test
RM_chi2 <- matrix(c(restric_stats_summary[system=="Restriction-modification" & presence==T]$nEffectif, 
                    absence =restric_stats_summary[system=="Restriction-modification" & presence==F]$nEffectif), ncol = 2 )
colnames(RM_chi2) <- c("present", "absence")
rownames(RM_chi2)  <- c("I", "II", "III")  
chisq.test(RM_chi2)
  
CRISPR_chi2 <- matrix(c(restric_stats_summary[system=="CRISPR" & presence==T]$nEffectif, 
                        absence =restric_stats_summary[system=="CRISPR" & presence==F]$nEffectif), ncol = 2 )
colnames(CRISPR_chi2) <- c("present", "absence")
rownames(CRISPR_chi2)  <- c("I", "II", "III")  
chisq.test(CRISPR_chi2)

