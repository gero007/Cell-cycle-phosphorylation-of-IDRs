# PSP Search Result Legend	 
# COLUMN HEADER	DESCRIPTION
# PROTEIN	Primary name in PhosphoSitePlusÂ® (PSP)
# GENE	Unique ID (UID) assigned to a gene by the HUGO Nomenclature Committee and adapted by HUPO for the cognate protein.
# ORGANISM	The species from which either the kinase and/or substrate proteins have been derived
# ACC_ID	Primary accession ID in PSP
# MW_(DA)	Molecular weight in Daltons.
# MOD_RSD	The location within the protein sequence of the amino acid residue that is posttranslationally modified (modsite).
# SITE_+/-7_AA	The modsite plus flanking sequence (+/- 7 rsds.). The modsite, as well as other residues within the flanking sequence that are known to be posttranslationally modified, are displayed in lower-case letters.
# SITE_GRP_ID	Unique identifier of the modification site and its homologous sites in all proteoforms and species
# DOMAIN	The Pfam residue in which the modsite is located.
# METHOD	The experimental method in which kinase-substrate relationships were reported. Vivo = in vivo: determined from reactions within living cells, cell cultures or organisms. Vitro = in vitro: determined from reactions outside of living cellular structures.
# LTP_LIT	The number of literature records using low-throughput (LTP) experimental techniques in which the modsite has been reported. LTP results may be more reliable than MS2 results.
# MS2_LIT	The number of literature records using tandem mass spectrometry (MS2) experimental techniques in which the modsite has been observed.
# MS2_CST	The number of curation sets (CS) derived from MS2 experiments performed at Cell Signaling Technology (CST) in which the modsite has been observed.
# CST_CAT#	The catalog number(s) of CST antibodies specific for the associated modsite.

#there's overlapping between the targets. I create 3 differents DF and I will merge them

phosphosite_CDK1Target_HS <- read_delim("PSP/PSP_CDK1_Target_HS.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_CDK2Target_HS <- read_delim("PSP/PSP_CDK2_Target_HS.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_CDK3Target_HS <- read_delim("PSP/PSP_CDK3_Target_HS.tab","\t", escape_double = FALSE, trim_ws = TRUE)

phosphosite_CDK1Target_HS <- phosphosite_CDK1Target_HS %>% rename(Kinase.CDK1 = KINASE)
phosphosite_CDK2Target_HS <- phosphosite_CDK2Target_HS %>% rename(Kinase.CDK2 = KINASE)
phosphosite_CDK3Target_HS <- phosphosite_CDK3Target_HS %>% rename(Kinase.CDK3 = KINASE)

phosphosite_CDKFamilyTarget_HS <- merge.data.frame(phosphosite_CDK1Target_HS,phosphosite_CDK2Target_HS,by = c("GENE","PROTEIN","ACC#","MOD_RSD"),all = T)
# colnames(phosphosite_CDKFamilyTarget_HS)[-(1:4)]<-c("METHOD.CDK1","LTP_LIT.CDK1","MS2_LIT.CDK1","MS2_CST.CDK1","Kinase.CDK1","METHOD.CDK2","LTP_LIT.CDK2","MS2_LIT.CDK2","MS2_CST.CDK2","Kinase.CDK2")
phosphosite_CDKFamilyTarget_HS <- merge.data.frame(phosphosite_CDKFamilyTarget_HS,phosphosite_CDK3Target_HS,by = c("GENE","PROTEIN","ACC#","MOD_RSD"),all = T)
# colnames(phosphosite_CDKFamilyTarget_HS)[-(1:14)]<-c("METHOD.CDK3","LTP_LIT.CDK3","MS2_LIT.CDK3","MS2_CST.CDK3","Kinase.CDK3")

# rows (11,12,13,14,15,16,49,79,80,151,167,168,170,190,206,211,320,321,322,323,324,327,358,412,439,440,443,444,445,446,447,448,449,507,508,509,510,531,532,533,534,596,633,762,834,849,850,851,882)
# having isoforms of the proteins. Convert them to the canonical form


###############################SOLVED############################

# vinexin iso10 = O60504_VAR_055019 rows(37)
# phosphosite S563
# No change in the sequence

phosphosite_CDKFamilyTarget_HS$`ACC#`[755]<- "O60504"
#------------------

# BIM iso2 = O43521-2 rows(62)
# change  42-101: Missing.
# phosphosite S44
phosphosite_CDKFamilyTarget_HS$`ACC#`[63]<- "O43521"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[63]<-"S104"

#------------------

# CDC25B iso2 = P30305-2 rows(117,118)
# change  68-81: Missing.
# phosphosite S146, S307

phosphosite_CDKFamilyTarget_HS$`ACC#`[125:126]<- "P30305"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[125]<-"S160"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[126]<-"S321"


#------------------

# CROCC iso2 = Q5TZA2-2 rows(594)
# change         1-697: Missing.     1984-1990: Missing.
# phosphosite S763

phosphosite_CDKFamilyTarget_HS$`ACC#`[192]<- "Q5TZA2"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[192]<-"S1460"

#------------------

# dUTPase iso2 = P33316-2 rows(214)
# change   1-93: MTPLCPRPAL...KAGGSPAPGP → MPCSE
# phosphosite S11

phosphosite_CDKFamilyTarget_HS$`ACC#`[224]<- "P33316"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[224]<-"S99"

#------------------

# EPB41 iso2 = P11171-2 rows(238,168)
# change   616-648: Missing.
# phosphosite T60, S679

phosphosite_CDKFamilyTarget_HS$`ACC#`[251:252]<- "P11171"
# T60 doesn't change
phosphosite_CDKFamilyTarget_HS$MOD_RSD[251]<-"S712"

#------------------

# FOXM1 iso2 = Q08050-2 rows(189)
# change    326-340: Missing.
# phosphosite T596

phosphosite_CDKFamilyTarget_HS$`ACC#`[281]<- "Q08050"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[281]<-"T611"

#------------------

# HMGA1 iso2 = P17096-2 rows(293)
# change    35-45: Missing.
# phosphosite T42

phosphosite_CDKFamilyTarget_HS$`ACC#`[310]<- "P17096"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[310]<-"T53"

#------------------

# ING1 iso2 = Q9UK53-2 rows(307)
# change    35-45: Missing.
# phosphosite S126

phosphosite_CDKFamilyTarget_HS$`ACC#`[324]<- "Q9UK53"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[324]<-"S269"

#------------------

# Tau iso8 = P10636-8 rows(376:382)
# change     125-375: Missing.     395-460: Missing.
# phosphosite  T153 S202 T205 T212 T231 S235 S404

phosphosite_CDKFamilyTarget_HS$`ACC#`[406:412]<- "P10636"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[409] <- "T470"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[406] <- "S519"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[410] <- "T522"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[411] <- "T529"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[412] <- "T548"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[407] <- "S552"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[408] <- "S721"

#------------------

# NUP98 iso2 = P52948-2 rows(493:497)
# change     393-409: Missing.     1502-1576: RHYDLNQLLE...FVLLHIDNSG → S
# phosphosite T529 T536 S595 S606 T653

phosphosite_CDKFamilyTarget_HS$`ACC#`[560:567]<- "P52948"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[562]<-"T546"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[563]<-"T553"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[560]<-"S612"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[561]<-"S623"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[564]<-"T670"

#------------------

# ODF2 iso3 = Q5BJF6-3 rows(500)
# change     1-41: MSASSSGGSP...PCGAPSVTVT → MKDRSSTPPL...LPKPSATSSQ.   65-83: Missing
# phosphosite S796

phosphosite_CDKFamilyTarget_HS$`ACC#`[567]<- "Q5BJF6"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[567]<-"S820"

#------------------

# PTPN2 iso2 = P17706-2 rows(554)
# change      382-415: WLYWQPILTKMGFMSVILVGAFVGWTLFFQQNAL → PRLTDT
# phosphosite S304

phosphosite_CDKFamilyTarget_HS$`ACC#`[624]<- "P17706"


#------------------


# AML1 iso8 = Q01196-8  rows(639:644)
#  change 1-5: MRIPV → MASDSIFESFPSYPQCFMRECILGMNPSRDVH, the difference is -27
# phosphosites "S48"   "S276"   "S293"   "T300"   "S303"   "S424"

phosphosite_CDKFamilyTarget_HS$`ACC#`[712:717]<- "Q01196"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[716]<-"S21"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[712]<-"S249"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[713]<-"S266"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[717]<-"T273"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[714]<-"S276"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[715]<-"S397"

#------------------

# SEPT9 iso2 = Q9UHD8-2 rows(650)
# change     1-25: MKKSYSGGTRTSSGRLRRLGDSSGP → MERDRIS
# phosphosite T24

phosphosite_CDKFamilyTarget_HS$`ACC#`[724]<- "Q9UHD8"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[724]<-"T38"

#------------------

# SIRT2 iso2 = Q8IXJ6-2 rows(663)
# change      1-37: Missing
# phosphosite S331

phosphosite_CDKFamilyTarget_HS$`ACC#`[737]<- "Q9UHD8"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[737]<-"S368"

#------------------

# SORBS3  iso10 = O60504_VAR_055019 rows(679)
# change      1canonical has an SP in that position. keeping the coordinates
# phosphosite S563
# 
# phosphosite_CDKFamilyTarget_HS$`ACC#`[679]<- "O60504"
# phosphosite_CDKFamilyTarget_HS$MOD_RSD[679]<-"S368"

#------------------

# SUN1  iso9 = O94901-9 rows(695:696)
# change      1canonical has an SP in that position. keeping the coordinates
# phosphosite S48 S333
#S48 doesn't change, S333 isn't present. removal below
phosphosite_CDKFamilyTarget_HS$`ACC#`[774:775]<- "O94901"

#------------------

# UHRF1  iso2 = Q96T88-2 rows(771)
# change      1-1: M → MGVFAVPPLSSDTM
# phosphosite S674

phosphosite_CDKFamilyTarget_HS$`ACC#`[852]<- "Q96T88"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[852]<-"S661"

#------------------

# VGLL4 iso4 = Q14135-4 rows(786:789)
# change      1-22: METPLDVLSRAASLVHADDEKR → MLFMKMDLLNYQYLDKMNNNIGILCYEG
# phosphosite  S58 S155 T159 S280

phosphosite_CDKFamilyTarget_HS$`ACC#`[867:870]<- "Q14135"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[869] <- "S52"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[867] <- "S149"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[870] <- "T153"
phosphosite_CDKFamilyTarget_HS$MOD_RSD[868] <- "S274"


#------------------

# ERK1 iso3 = P27361-3 rows(170)
# change  340-379: PVAEEPFTFAMELDDLPKERLKELIFQETARFQPGVLEAP → VGQSPAAVGLGAGEQGGT
# phosphosite S343
# The paper states that phosphorylation only happens in the isoform 3, even if they have evidence of phosphorylation in other sites. They don't investigate other posibles phosphosites
# removal below


phosphosite_CDKFamilyTarget_HS<- phosphosite_CDKFamilyTarget_HS[c(-401,-774),]



#################################################################

# Remmove the 'iso' part of the names
phosphosite_CDKFamilyTarget_HS$PROTEIN <- str_split_fixed(phosphosite_CDKFamilyTarget_HS$PROTEIN," iso",n = 2)[,1]

human_data <- phosphosite_CDKFamilyTarget_HS %>% unite(col = "KINASE", Kinase.CDK1,Kinase.CDK2,Kinase.CDK3,sep = ",",remove = T,na.rm=T )
# human_data <- human_data %>% unite(col = "method", METHOD.CDK1,METHOD.CDK2,METHOD.CDK3,sep = ",",remove = T,na.rm=T )


human_data <- human_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at(c("MOD_RSD","KINASE"),function(x){paste(x, collapse=",")})

human_data$MOD_RSD <- sapply(human_data$MOD_RSD, function(x){ return(paste(substr(strsplit(x,",")[[1]],start = 2,stop = 2000),collapse = ","))})
human_data$KINASE <- sapply(human_data$KINASE, function(x){ return(paste(unique(strsplit(x,",")[[1]]),collapse = ","))})
# algunos valores tienen espacios, por eso el gsub
# human_data$method <- sapply(human_data$method, function(x){ return(paste(unique(gsub(" ", "", strsplit(x,",")[[1]], fixed = TRUE)),collapse = ","))})



# The data is exported and mannually curated afterwards. Psites from bibliography were added. All the information is then deposited in the file human_data_curated.tab