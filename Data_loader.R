library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggsci)
library(scales)
library(reshape2)

source("functions/expectedVsObserved.R")

######IUPred data loading########
disoPath <- "predictions/IUpred_run_human/"
file.names <- dir(disoPath, pattern =".iupred")
ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
IUPred_predictions <- data.frame(ids)
colnames(IUPred_predictions) <- "ID"
IUPred_predictions$ID <- as.character(IUPred_predictions$ID)
# positions <- list()
scores <- list()
IUPred_disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","IUPred_SCORE")
  aux_table$IUPred_DISO <- aux_table$IUPred_SCORE>=0.5
  IUPred_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  # positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$IUPred_SCORE)
  # IUPred_predictions[n,"phospho"] <- NULL
  # IUPred_predictions[n,"threshold"] <- 0.5
  IUPred_disordered[[n]] <- which(aux_table$IUPred_DISO)
}
# IUPred_predictions$positions <-positions
IUPred_predictions$IUPredScores <-scores
IUPred_predictions$IUPred_disordered <-IUPred_disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% IUPred_predictions$ID]
IUPred_predictions$length <- nchar(IUPred_predictions$sequence)

########SPOT###########
disoPath <- "predictions/spot_run_human/"
file.names <- dir(disoPath, pattern =".spotd")
ids <-unlist(strsplit(x = file.names,split = ".spotd",fixed = T))
SPOT_predictions <- data.frame(ids)
colnames(SPOT_predictions) <- "ID"
SPOT_predictions$ID <- as.character(SPOT_predictions$ID)
# positions <- list()
scores <- list()
SPOT_disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","SPOT_SCORE","SPOT_DISO")
  SPOT_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  # positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$SPOT_SCORE)
  # IUPred_predictions[n,"phospho"] <- NULL
  # IUPred_predictions[n,"threshold"] <- 0.5
  SPOT_disordered[[n]] <- which(aux_table$SPOT_DISO=="D")
}
# SPOT_predictions$positions <-positions
SPOT_predictions$SPOT_scores <-scores
SPOT_predictions$SPOT_disordered <-SPOT_disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% IUPred_predictions$ID]
SPOT_predictions$length <- nchar(SPOT_predictions$sequence)


#Merge disorder predictions

disorder_predictions <- merge.data.frame(IUPred_predictions,SPOT_predictions,by = c("ID","sequence","length"))

# clean the WS
rm(IUPred_disordered,SPOT_disordered,scores)



human_data <- read_delim("PSP/human_data_curated_V2_cleaned.tab", 
                         "\t", escape_double = FALSE, col_types = cols(MOD_RSD = col_character()), 
                         trim_ws = TRUE)

human_data$psites_CDK1 <- lapply(human_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
human_data$target <- rep("Cdk1 target",nrow(human_data))
human_data$target <- as.factor(human_data$target)

human_universe_data <- read_delim("PSP/Phosphorylation_site_dataset", 
                                  "\t", escape_double = FALSE, col_types = cols(HU_CHR_LOC = col_skip(), 
                                                                                SITE_GRP_ID = col_skip(), MW_kD = col_skip(), 
                                                                                DOMAIN = col_skip(), `SITE_+/-7_AA` = col_skip(), 
                                                                                LT_LIT = col_skip(), MS_LIT = col_skip(), 
                                                                                MS_CST = col_skip(), `CST_CAT#` = col_skip()), 
                                  trim_ws = TRUE)

human_universe_data<-rename(human_universe_data,c(`ACC#`=ACC_ID))
# Select human proteins
human_universe_data <- subset(human_universe_data, ORGANISM == "human")
# Remove all information related to isoforms
human_universe_data <- subset(human_universe_data, !grepl("-",`ACC#`))
human_universe_data <- subset(human_universe_data, !grepl(" iso[0-9]",`PROTEIN`))
# Select targets with pS or pT
# human_universe_data<-subset(human_universe_data,(substr(MOD_RSD,1,1)=="S"|substr(MOD_RSD,1,1)=="T"))
# format the MOD_RSD column and group by gene/protein/uniprot
human_universe_data <- human_universe_data %>% mutate(MOD_RSD=substr(MOD_RSD,1,nchar(MOD_RSD)-2)) 
human_universe_data <- human_universe_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at("MOD_RSD",function(x){paste(substr(x,2,2000), collapse=",")})
# Generate the psite column, with the vectors containing psites (for the contingency table analysis)
human_universe_data <- as.data.frame(human_universe_data)
human_universe_data$psites <- lapply(human_universe_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})

# merge the tables and mark CDK1 targets and non CDK1 targets. Only by ACC, protein names in human data are still with  the isoform nomenclature
human_data <- merge.data.frame(x = human_data,y = human_universe_data,by = c("ACC#"),all = T,suffixes = c("_CDK1","_ALL"))
#remove PROTEIN and GENE comumns from x, and rename the y columns
human_data$GENE_CDK1<-NULL
human_data$PROTEIN_CDK1<-NULL
human_data <- rename(human_data,c(PROTEIN=PROTEIN_ALL,GENE=GENE_ALL))
# Adding the target category "Non Cdk1 target"
levels(human_data$target) <- c("Cdk1 target","Non Cdk1 target")
human_data$target[is.na(human_data$target)]<-"Non Cdk1 target"

human_data<-merge.data.frame(human_data,disorder_predictions,by.x = "ACC#",by.y = "ID")

human_data$psites_count <- sapply(human_data$psites, length)
human_data$psites_CDK1_count <- sapply(human_data$psites_CDK1, length)

human_data$ST_residues <- lapply(gregexpr("S|T",(human_data$sequence)), FUN=as.numeric)


#########IUPred#########
IUPred_expVSobs <- expectedVsObserved(human_data,"IUPred")

human_data$IUPred_psites_obsv_diso <- IUPred_expVSobs$phosphoDiso_obs
human_data$IUPred_psites_expct_diso <- IUPred_expVSobs$phosphoDiso_expct
human_data$IUPred_psites_expct_diso_prob <- IUPred_expVSobs$phosphoDiso_expct_prob


human_data <- as.data.table(human_data)
human_data[,IUPred_binom := purrr::pmap(.(IUPred_psites_obsv_diso, psites_count, IUPred_psites_expct_diso_prob), binom.test, alternative="greater")]
human_data[,IUPred_binom_p := mapply("[[", IUPred_binom, "p.value", SIMPLIFY = T)]
human_data[,IUPred_binom_q := mapply(p.adjust, IUPred_binom_p)]
human_data[,IUPred_binom_sig := factor(ifelse(IUPred_binom_q < 0.05, ifelse(IUPred_binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]


#########SPOT#########
SPOT_expVSobs <- expectedVsObserved(human_data,"SPOT")

human_data$SPOT_psites_obsv_diso <- SPOT_expVSobs$phosphoDiso_obs
human_data$SPOT_psites_expct_diso <- SPOT_expVSobs$phosphoDiso_expct
human_data$SPOT_psites_expct_diso_prob <- SPOT_expVSobs$phosphoDiso_expct_prob


human_data <- as.data.table(human_data)
human_data[,SPOT_binom := purrr::pmap(.(SPOT_psites_obsv_diso, psites_count, SPOT_psites_expct_diso_prob), binom.test, alternative="greater")]
human_data[,SPOT_binom_p := mapply("[[", SPOT_binom, "p.value", SIMPLIFY = T)]
human_data[,SPOT_binom_q := mapply(p.adjust, SPOT_binom_p)]
human_data[,SPOT_binom_sig := factor(ifelse(SPOT_binom_q < 0.05, ifelse(SPOT_binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]

save.image("plotting_data_test.RData")
