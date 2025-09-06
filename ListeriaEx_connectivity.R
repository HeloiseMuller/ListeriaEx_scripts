library(data.table)
library(dplyr)
library(igraph)

setwd("project_Listeria/")

#Choose the dataset to apply the script on
datatset = "global"
#datatset = "local"

if(dataset=="global"){
  allMGE <- fread("allComplete_AND_allhybrids_MGE_filtered_nodes.tbl")
  net <- read_graph("allComplete_AND_allhybrids_MGE_filtered.gml", format = "gml")
  membershipsGiant = c(1,2,5,6,7)
  corresMembership = data.table(inNetwork=membershipsGiant, inFigure = c("B", "C", "1", "3", "2"))
} 

if(datatset=="local"){
  allMGE <- fread("all_936SRR_MGE_filtered_nodes.tbl")
  allMGE <- rename(allMGE, MGEname="nodeName")
  net <- read_graph("graph_all_936SRR_MGE_filtered.gml", format = "gml")
  membershipsGiant = c(1,9,10)
  corresMembership = data.table(inNetwork=membershipsGiant,  inFigure = c("2", "1", "3"))
  
}


###########################

#Rest of the code is the same regardless of the dataset
dt <- data.table(as_data_frame(net))
dt <- rename(dt, QID="from", SID="to")

dt <- left_join(dt, allMGE[, c("MGEname", "membership", "type")], by = c("QID"="MGEname")) %>%  left_join(., allMGE[, c("MGEname", "membership", "type")], by = c("SID"="MGEname")) 

#Keep only connections inside the giant component
dt_giant <- dt[membership.x %in% membershipsGiant & membership.y %in% membershipsGiant,]

#Write numbering as in figure
dt_giant <- left_join(dt_giant, corresMembership, by = c("membership.x"="inNetwork")) %>% left_join(., corresMembership, by = c("membership.y"="inNetwork"))

dt_giant[, sameCom := ifelse(inFigure.x==inFigure.y, T, F)]


#Count number of connections per nodes per types
halfOne <- dt_giant[, c( "QID", "inFigure.x", "inFigure.y", "sameCom", "type.x", "type.y")] 
halfTwo <-   dt_giant[, c( "SID", "inFigure.y", "inFigure.x", "sameCom", "type.y", "type.x")]
names(halfOne) = names(halfTwo) =  c("ID", "community", "communityCo", "sameCom", "type", "typeCo")
dt_giant_perCO <- rbind(halfOne, halfTwo)

coPerNode <- dt_giant_perCO%>% group_by(ID, type) %>% summarize(nTot = n()) %>% data.table

coPerPaired <- dt_giant_perCO %>% group_by(ID, type, typeCo) %>% summarize(nTot = n()) %>% data.table %>%
  mutate(., pairType = ifelse(type<typeCo, paste0(type, "-",typeCo), paste0(typeCo, "-", type)))

#Also need to add nodes with 0 connections
nodesGiant <- allMGE[membership %in% membershipsGiant, c("MGEname", "membership", "type")]

#I give 4 lines per node: a connection with transposon, one with phage, one with plasmid, one with PP
dummy <- data.table(ID = nodesGiant$MGEname, type =  nodesGiant$type, typeCo = "transposon") %>% 
  rbind(., data.table(ID = nodesGiant$MGEname, type =  nodesGiant$type, typeCo = "phage")) %>%
  rbind(., data.table(ID = nodesGiant$MGEname, type =  nodesGiant$type, typeCo = "plasmid")) %>%
  rbind(., data.table(ID = nodesGiant$MGEname, type =  nodesGiant$type, typeCo = "PP")) %>%
  mutate(., pairType = ifelse(type<typeCo, paste0(type, "-",typeCo), paste0(typeCo, "-", type))) %>%
  filter(., pairType != "phage-transposon")

#Add co that don't exist
coPerPaired <- left_join(dummy, coPerPaired)
coPerPaired[is.na(nTot),]$nTot = 0

coPerPaired <- mutate(coPerPaired, nType = ifelse(type=="phage", 2171, ifelse(type=="transposon", 89, 170)), 
                      nTypeCo = ifelse(typeCo =="phage", 2171, ifelse(typeCo=="transposon", 89, 170))) %>%
  mutate(., nCo_relative2 = nTot / (nType-1))

### stat per com
coPerPaired_com <- dt_giant_perCO %>% group_by(ID, type, typeCo, community, communityCo, sameCom) %>% summarize(nTot = n()) %>% data.table %>%
  mutate(., pairType = ifelse(type<typeCo, paste0(type, "-",typeCo), paste0(typeCo, "-", type)))

coPerPaired_com[, com_type_Co:=paste(communityCo, typeCo, sep = ":")]

#Count number of nodes per com
nNodes_perCom_perType <- dt_giant_perCO[, c("ID", "community", "type")] %>% unique %>% group_by(community, type) %>% summarise(nNode = n()) %>% data.table
nNodes_perCom_perType[, com_type:=paste(community, type, sep = ":")]

#Max number of nodes with same type and same com:
nNodes_perCom_perType[, maxCo_sameCom_sameType := nNode-1]


#Calculate proportion out of possibilities
coPerPaired_com <-
  rbind(
    rename(left_join(coPerPaired_com[sameCom==T & type==typeCo,],
                     nNodes_perCom_perType[, c("com_type", "maxCo_sameCom_sameType")],
                     by = c("com_type_Co" = "com_type")), maxCo="maxCo_sameCom_sameType"),
    rename(left_join(coPerPaired_com[sameCom!=T | type!=typeCo,],
                     nNodes_perCom_perType[, c("com_type", "nNode")],
                     by = c("com_type_Co" = "com_type")), maxCo="nNode")
  )

coPerPaired_com$proportion <- coPerPaired_com$nTot/coPerPaired_com$maxCo

#Will show tot nuber of node of a com in panel
nNode_legend_panel <- distinct(dt_giant_perCO[, c("ID", "community")]) %>% group_by(community) %>% summarize(nTot_com=n()) %>% data.table()

#OR can also show more specificaly nymber of node of each type
nNode_legend_panel2 <- distinct(dt_giant_perCO[, c("ID", "community", "type")]) %>% group_by(community, type) %>% summarize(nTot_com_type=n()) %>% data.table()
setorder(nNode_legend_panel2, -nTot_com_type)
nNode_legend_panel2 <- mutate(nNode_legend_panel2, type_short = ifelse(type=="phage", "Ph", ifelse(type=="transposon", "Tn", ifelse(type=="plasmid", "Pl",ifelse(type=="PP", "PP",NA))))) %>%
  mutate(nNode_legend_panel2, nTot_com_type = paste0(nTot_com_type, " ", type_short) ) %>% group_by(community) %>% summarize(nTot_com_types = paste(nTot_com_type, collapse = " + "))


coPerPaired_com <- left_join(coPerPaired_com, nNode_legend_panel) %>% mutate(., panel_title = paste0("Community ", community, " (n = ", nTot_com, ")")) %>% select(., -nTot_com)
coPerPaired_com <- left_join(coPerPaired_com, nNode_legend_panel2) %>% mutate(., panel_title = paste0("Community ", community, " (n = ", nTot_com_types, ")")) %>% select(., -nTot_com_types)

#MAke the figure
p <- ggplot(coPerPaired_com, aes(x=paste0("Community ", communityCo), y=proportion)) +
  geom_boxplot(fill="white", outlier.colour=NA, position=position_dodge(width=0.9), aes(col=pairType)) +
  geom_point(alpha = .5, position = position_dodge(width=0.9), aes(col = pairType)) + 
  #facet_wrap(~paste0("Community ", community)) +
  facet_wrap(~panel_title) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  ylab("Proportion of connections with the MGEs of each community") +
  scale_color_manual("Connection type",
                     values = c("phage-phage" = "#f6dfd1ff",
                                "phage-plasmid"="#d9689eff", 
                                "phage-PP"="#252525",
                                "plasmid-plasmid"="#cbcbcbff",
                                "plasmid-PP"="#252525",
                                "plasmid-transposon" = "#8fc379ff",
                                "PP-PP" = "#FEB24C", 
                                "PP-transposon"="#8C2D04",
                                "transposon-transposon" = "#6292c3ff")
  ) +
  theme(legend.position=c(.85,.15))

#I should add how many nodes are ncluded in each categories:
nNode_legend <- coPerPaired_com %>% group_by(community, communityCo, pairType) %>% summarize(n=n()) %>% data.table
#Or say more specificaly what MGE it is:
nNode_legend <- coPerPaired_com %>% group_by(community, communityCo, type, pairType) %>% summarize(n=n()) %>% data.table

#Add values of 0
a <- distinct(nNode_legend[, c("communityCo", "pairType")])
b <- cbind(a, community = corresMembership[1,]$inFigure) %>% 
  rbind(., cbind(a, community = corresMembership[2,]$inFigure)) %>% 
  rbind(., cbind(a, community = corresMembership[3,]$inFigure)) %>% 
  rbind(., cbind(a, community = corresMembership[4,]$inFigure)) %>% 
  rbind(., cbind(a, community = corresMembership[5,]$inFigure)) 
nNode_legend <- left_join(b, nNode_legend)
nNode_legend[is.na(n),]$n = 0


#Do not show intra com since annoying when several types of MGE
nNode_legend <- nNode_legend[community!=communityCo]

#Avoid a value to overlap 0. In certain case these lines could be a problem,
#but since i  never have T connection type with an inter com, I am good
setorder(nNode_legend, -n)
nNode_legend <- nNode_legend[!duplicated(paste(nNode_legend$community, nNode_legend$communityCo))]

#Since i now use panel_title, I need the same thing here
nNode_legend <- left_join(nNode_legend, nNode_legend_panel) %>%  mutate(., panel_title = paste0("Community ", community, " (n = ", nTot_com, ")")) %>% select(., -nTot_com)
#Or 2nd version:
nNode_legend <- left_join(nNode_legend, nNode_legend_panel2) %>%  mutate(., panel_title = paste0("Community ", community, " (n = ", nTot_com_types, ")")) %>% select(., -nTot_com_types)

nNode_legend <- mutate(nNode_legend, type_short = ifelse(type=="phage", "Ph", ifelse(type=="transposon", "Tn", ifelse(type=="plasmid", "Pl", ifelse(type=="PP", "PP", NA))))) %>%
  mutate(., n = ifelse(is.na(type), n, paste0(n, " ", type_short)))

#######################
# Save the figure

##Figure local dataset
if(dataset=="local"){  
  #pdf("FigureS3.pdf", width = 10, height = 6)
  p + geom_text(data =  nNode_legend,aes(y=-0.1, label = n), size = 3, position = position_dodge(width = 1))  
  dev.off()
}

##Figure global dataset
if(dataset=="global"){
  pdf("FigureS6.pdf", width = 10, height = 6)
  p + geom_text(data =  nNode_legend, aes(y=-0.1, label = n), size = 3, position = position_dodge(width = 1))  
  dev.off()
}