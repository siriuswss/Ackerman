### R code from vignette source 'analysis_Ackerman_26April2013.Rnw'

###################################################
### code chunk number 1: analysis_Ackerman_26April2013.Rnw:24-29
###################################################
  library(stringr)
  UserID = Sys.getenv("LOGNAME")
  WorkDir = paste(getwd(),"/",sep="")
  Rversion = paste(str_trim(substr(sessionInfo()[1]$R.version$version.string,1,17), side="right"),","," ",UserID,sep="",collapse=NULL)
  len = (nchar(Rversion)/100) + ((nchar(Rversion)/100)/4)


###################################################
### code chunk number 2: load-packages
###################################################
rm(list=ls(all=TRUE))
setwd("N:/cavd/VIMCs/Antibody/FcR Comparison/analysis/FinalStudy/code")
#Source("../macro/load_libraries.R")
#load_libraries()
library(ggplot2)


###################################################
### code chunk number 3: ReadIn_PrepareData
###################################################
code = read.csv("../encryption_Ackerman/encryption.csv")
colnames(code) = tolower(colnames(code))


###################################################
### code chunk number 4: ReadIn_PrepareData1_Ack
###################################################
man = read.csv("../qdata/FcRFinal_Ackerman_20130314.csv")
colnames(man) = tolower(colnames(man))
man = man[,c('assay_subtype','detectionreagent','bead_no','antigen_name','well_role','spec_type','sample','study_control_name',
             'name','standard_name','expconc','starting_concentration','dilution_factor','fi','fibkgd')]

man$assay_subtype = interaction("Acker", man$assay_subtype)

# Get all obs where well_role is Sample
man = man[man$well_role=="Unknown",]

# Put name into sample
man$sample = as.character(man$sample)
man$name = as.character(man$name)
man$sample = ifelse(man$name!="",man$name,man$sample)

# Sort data
man = man[with(man, order(sample,antigen_name,detectionreagent)),]

# Remove pilot study samples FCAC001 to FCAC010
man = man[!man$sample %in% c("FCAC001","FCAC002","FCAC003","FCAC004","FCAC005","FCAC006","FCAC007","FCAC008","FCAC009","FCAC010"),]

# Merge coding
man = merge(man, code, by.x=c("sample"),by.y="blinded_sample_id",all.x=TRUE)
man$group_assignment = as.character(man$group_assignment)
man$study_control_name = as.character(man$study_control_name)
man$group_assignment = ifelse(man$study_control_name=="IVIG", "IVIG", man$group_assignment)
man$group_assignment = ifelse(man$study_control_name=="HIVIG", "HIVIG & HIVIG-C", man$group_assignment)
man$group_assignment = ifelse(man$study_control_name=="HIVIG-C", "HIVIG & HIVIG-C", man$group_assignment)

# Change colnames
names(man)[names(man)=="assay_subtype"] = "Assay"
names(man)[names(man)=="group_assignment"] = "Group"

# Remove plasma samples
out = data.frame( do.call( rbind, strsplit( man$sample, ' ' ) ) ) 
names(out) = paste('column',1:2,sep="")
man = cbind(man, out)
man = man[man$column1!="Plasma",][,c(1:19)]

# Remove rows where Group is Blank
man = man[man$Group!="",]
# Order levels for Group variable
man$Group <- factor(man$Group, levels = c('CHRONIC Treated','chronic untreated','elite controller','VAX004','Placebo','HIV Seroneg','IVIG','HIVIG & HIVIG-C'),ordered = TRUE)

rm(out)


###################################################
### code chunk number 5: ReadIn_MetaData_Ack
###################################################
meta = read.csv("../qdata/FcRPilotStudy_AntigenMetadata_20130320.csv")
colnames(meta) = tolower(colnames(meta))
meta = meta[,c('antigenname_ackermanlab','beadnum_ackerman_fcr_final','study.specific.antigen.name','subtype')]

meta$aginfo = interaction(meta$study.specific.antigen.name, meta$subtype)
meta$aginfo = as.character(meta$aginfo)

# Conflict in HSV-1 and HSV-2. So resolve by making it HSV. 
meta$aginfo = ifelse(meta$aginfo=="HSV-2 gG.", "HSV gG.", meta$aginfo)

# Missing Study_Specific_Antigen_Name for "gp120 PVO"
meta$aginfo = ifelse(meta$aginfo==".B", "r.gp120 PVO.B", meta$aginfo)

# Merge man and meta
man = merge(man, meta, by.x=c("bead_no"),by.y="beadnum_ackerman_fcr_final",all.x=TRUE)


###################################################
### code chunk number 6: analysis_Ackerman_26April2013.Rnw:249-250
###################################################
xtable(cbind(c("Total IgA", "Total IgG", rep("FcR",5),rep("IgG Subclass",4),rep("Lectin",5)),as.character(sort(unique(man$detectionreagent)))))


###################################################
### code chunk number 7: analysis_Ackerman_26April2013.Rnw:253-259
###################################################
aglist_gp120 <- grep("gp120",unique(man$antigen_name),ignore.case=T,value=T)
aglist_gp140 <- grep("gp140",unique(man$antigen_name),ignore.case=T,value=T)
aglist_gp41 <- grep("gp41",unique(man$antigen_name),ignore.case=T,value=T)
aglist_SOSIP <- grep("SOSIP",unique(man$antigen_name),ignore.case=T,value=T)

xtable(cbind(c("Total IgA", "Total IgG", rep("FcR",5),rep("IgG Subclass",4),rep("Lectin",5)),as.character(sort(unique(man$detectionreagent)))))


###################################################
### code chunk number 8: BoxPlot1_1_Ack
###################################################

data = ahuigg = man[man$detectionreagent=="ahuIgG",]

q = ggplot(data, aes(x=factor(Group), y=fi)) 
q + theme_bw()  + ylab("FI") + xlab("Groups") +
  #scale_y_log10(breaks=c(0.05, 1,10,50), labels=c(eval(cutoffexp),1,10,50)) +
  scale_y_log10() +
  #geom_point(size=2,position = position_jitter(width = 0.2, height=0),aes(colour = factor(response)))  +
  geom_point(size=2,position = position_jitter(width = 0.2, height=0), colour="#FF9900")  +
  geom_boxplot(outlier.colour = "NA",alpha=0) +  
  opts(plot.margin = unit(c(1,1,0,1), "cm"),legend.position="bottom",
       legend.title=theme_text(size=7,vjust=-1.6),
       legend.text=theme_text(size=7,vjust=-1.6),
       axis.text.x = theme_text(size=9,angle = 45, hjust = 1,vjust=0.95), 
       axis.text.y = theme_text(size=8), 
       axis.title.x = theme_text(size=12), 
       axis.title.y = theme_text(angle=90, size=12,vjust=-0.4), 
       strip.text.x = theme_text(size = 7), 
       strip.text.y = theme_text(size = 12, face="bold", angle=270))      



###################################################
### code chunk number 9: BoxPlot1_2_Ack
###################################################

q = ggplot(data, aes(x=factor(Group), y=fibkgd)) 
q + theme_bw()  + ylab("FIBkgd") + xlab("Groups") +
  #scale_y_log10(breaks=c(0.05, 1,10,50), labels=c(eval(cutoffexp),1,10,50)) +
  scale_y_log10() +
  #geom_point(size=2,position = position_jitter(width = 0.2, height=0),aes(colour = factor(response)))  +
  geom_point(size=2,position = position_jitter(width = 0.2, height=0), colour="#FF9900")  +
  geom_boxplot(outlier.colour = "NA",alpha=0) +  
  opts(plot.margin = unit(c(1,1,0,1), "cm"),legend.position="bottom",
       legend.title=theme_text(size=7,vjust=-1.6),
       legend.text=theme_text(size=7,vjust=-1.6),
       axis.text.x = theme_text(size=9,angle = 45, hjust = 1,vjust=0.95), 
       axis.text.y = theme_text(size=8), 
       axis.title.x = theme_text(size=12), 
       axis.title.y = theme_text(angle=90, size=12,vjust=-0.4), 
       strip.text.x = theme_text(size = 7), 
       strip.text.y = theme_text(size = 12, face="bold", angle=270))

rm(data)


###################################################
### code chunk number 10: CrossAssay1_3_Ack
###################################################

#pdf(file="../graph/corr_heatmap120.pdf", height=12, width=12)
par(oma=c(10,0,0,10))
par(mar=c(0,0,0,0))
#dev.off()


