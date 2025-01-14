library(rtracklayer)
library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tximport)
library(clusterProfiler)
library(DESeq2)
library(limma)
library(viridis)
library(ggpubr)
library(edgeR)
library(ggrepel)
library(DOSE)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","ruthasinger_github","Opto-CLIP_v2","Figure6_Figure7_FigureS2_FigureS3")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#reading in GTF for transcript information
mm10_gtf=rtracklayer::import(file.path(Datadirectory,"mm10.ensGene.gtf.gz"))
mm10_gtf_df=mm10_gtf %>% as_tibble
mm10_Tx=mm10_gtf_df %>% dplyr::filter(type=="transcript")
mm10_Tx_final=mm10_Tx %>% dplyr::select(seqnames,start,end,width,strand,gene_id,transcript_id)
mm10_Tx_final$gene_name<- mapIds(org.Mm.eg.db,
                                 keys = mm10_Tx_final$gene_id,
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
mm10_Tx_final$entrez<- mapIds(org.Mm.eg.db,
                              keys = mm10_Tx_final$gene_id,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

mm10_Tx_final_noNAs=subset(mm10_Tx_final, is.na(gene_name) == FALSE)

################################################################################
#Bring in Opto-RiboTag data
RiboTag_Metadata <- read.delim(file.path(Datadirectory,"RiboTag_Metadata_justIP.txt"))
RiboTag_files <- file.path(Datadirectory,RiboTag_Metadata$Filename, "quant.sf")
names(RiboTag_files) <- RiboTag_Metadata$Filename
all(file.exists(RiboTag_files))

RiboTag_txi.tx <- tximport(RiboTag_files, type = "salmon", txOut = TRUE)
rownames(RiboTag_Metadata) <- RiboTag_Metadata$Filename

RiboTag_salmon_TPM=as.data.frame(RiboTag_txi.tx$abundance)
RiboTag_salmon_TPM <- tibble::rownames_to_column(RiboTag_salmon_TPM, "trancript_ID")
colnames(RiboTag_salmon_TPM) <- c("transcript_ID",paste(colnames(RiboTag_salmon_TPM[2:ncol(RiboTag_salmon_TPM)]), "TPM", sep = "."))

#output file of TPM counts for GEO
write.csv(RiboTag_salmon_TPM, file = file.path(Outdirectory,"RiboTag_salmon_TPM_forGEO.csv"), row.names=FALSE)

#Deseq
RiboTag_Metadata <- mutate_at(RiboTag_Metadata, vars(Group), as.factor)

RT_dds <- DESeqDataSetFromTximport(RiboTag_txi.tx, RiboTag_Metadata, ~ Group)

keep <- rowSums(counts(RT_dds)) >= 1
RT_dds <- RT_dds[keep, ]
RT_dds <- DESeq(RT_dds)

res=results(RT_dds, contrast=c("Group", "ChR2","Control"))
res=res[order(res$padj), ]
resdata=as.data.frame(res)
resdata <- tibble::rownames_to_column(resdata, "Transcript_ID")

RiboTag_resdata_TPM=merge(resdata,RiboTag_salmon_TPM, by.x="Transcript_ID", by.y="transcript_ID", all.x=TRUE, all.y=FALSE)

#add in gtf info
RiboTag_mm10_gtf_merged <- merge(RiboTag_resdata_TPM, mm10_Tx_final_noNAs, by.x="Transcript_ID", by.y="transcript_id")

#should be true
length(unique(RiboTag_mm10_gtf_merged$Transcript_ID)) == nrow(RiboTag_mm10_gtf_merged)

#should be false
length(unique(RiboTag_mm10_gtf_merged$gene_name)) == nrow(RiboTag_mm10_gtf_merged)

nrow(RiboTag_mm10_gtf_merged) #60270

RiboTag_mm10_gtf_merged_exp=RiboTag_mm10_gtf_merged %>% 
  dplyr::filter(RiboTag_ChR2_30min_IP_rep1_salmon.TPM > 0 & 
                RiboTag_ChR2_30min_IP_rep2_salmon.TPM > 0 & 
                RiboTag_ChR2_30min_IP_rep3_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep1_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep2_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep3_salmon.TPM > 0)

nrow(RiboTag_mm10_gtf_merged_exp) #30157

RiboTag_mm10_gtf_merged_exp= RiboTag_mm10_gtf_merged_exp %>%
  rowwise() %>%
  mutate(avg_rest_RT_TPM    = mean(c_across(c("RiboTag_Control_30min_IP_rep1_salmon.TPM", "RiboTag_Control_30min_IP_rep2_salmon.TPM","RiboTag_Control_30min_IP_rep3_salmon.TPM")), na.rm=TRUE),
         avg_opto_RT_TPM    = mean(c_across(c("RiboTag_ChR2_30min_IP_rep1_salmon.TPM", "RiboTag_ChR2_30min_IP_rep2_salmon.TPM","RiboTag_ChR2_30min_IP_rep3_salmon.TPM")), na.rm=TRUE),
         avg_both_RT_TPM    = mean(c_across(c("RiboTag_Control_30min_IP_rep1_salmon.TPM","RiboTag_Control_30min_IP_rep2_salmon.TPM","RiboTag_Control_30min_IP_rep3_salmon.TPM","RiboTag_ChR2_30min_IP_rep1_salmon.TPM", "RiboTag_ChR2_30min_IP_rep2_salmon.TPM","RiboTag_ChR2_30min_IP_rep3_salmon.TPM")), na.rm=TRUE))

#subset RiboTag data so it is a single Tx per gene- doing highest average across all samples
RiboTag_mm10_gtf_merged_highestExp=RiboTag_mm10_gtf_merged_exp %>%
  dplyr::group_by(gene_name) %>% 
  dplyr::arrange(desc(avg_both_RT_TPM)) %>%
  dplyr::filter(row_number()==1) 

nrow(RiboTag_mm10_gtf_merged_highestExp) #14353

#now should be true
length(unique(RiboTag_mm10_gtf_merged_highestExp$gene_name)) == nrow(RiboTag_mm10_gtf_merged_highestExp)

################################################################################################################################################################
#Bring in FMRP CLIP data
FMRP_CLIP_STAR_files=grep("transcriptome.unique.bed$",list.files(Datadirectory),value=TRUE)

FMRP_CLIP_STAR_files

transcripts=NULL
for (i in 1:(length(FMRP_CLIP_STAR_files))) {
  transcripts[[i]]=read_delim(file.path(Datadirectory,FMRP_CLIP_STAR_files[[i]]), delim = "\t", col_names = FALSE)
  colnames(transcripts[[i]])= c("Transcript_id", "Txstart", "Txend", "read_ID", "Txscore", "Txstrand")
  names(transcripts)[[i]]=vapply(strsplit(FMRP_CLIP_STAR_files[[i]], ".", fixed = TRUE), "[", "", 1)
}

All.FMRP.tags=do.call("rbind", transcripts)

All.FMRP.tags$sourcefile=rownames(All.FMRP.tags)
row.names(All.FMRP.tags) <- NULL

All.FMRP.tags= All.FMRP.tags %>% separate(sourcefile, into = c("Sample","two"), sep = "\\.", remove = TRUE) 
All.FMRP.tags= All.FMRP.tags %>% dplyr::select(!two)
nrow(All.FMRP.tags) #818547

#subset FMRP tags for those found in RiboTag transcripts
counts <- All.FMRP.tags %>%
  dplyr::filter(Transcript_id %in% RiboTag_mm10_gtf_merged_highestExp$Transcript_ID) %>%
  group_by(Transcript_id, Sample) %>%
  summarise(tags = dplyr::n()) %>% ungroup()

sum(counts$tags) #385172

counts2=counts %>% 
  separate(Sample, into = c("RBP","Treatment","rep"), sep = "\\_", remove = FALSE) %>%
  group_by(Transcript_id, Treatment) %>%                                                            
  summarise(BC = n())  %>%
  pivot_wider(names_from = Treatment, values_from = BC, names_prefix="BC_") %>% 
  replace(is.na(.), 0)

counts_wide=counts %>% 
  pivot_wider(names_from= Sample, values_from= tags) %>%
  replace(is.na(.), 0)

#output file of FMRP CLIP counts for GEO
write.csv(counts_wide, file = file.path(Outdirectory,"FMRP_counts_forGEO.csv"), row.names=FALSE)

counts_BC <- cbind(counts_wide, counts2[,2:3])

#filtering for BC
counts_BC4=counts_BC %>% filter(BC_ChR2 > 3 | BC_Control > 3)

nrow(counts_BC4) #4552

counts_BC4=counts_BC4 %>% 
  rowwise() %>%
  mutate(ChR2_max_tags=max(c_across(c("FMRP_ChR2_rep1","FMRP_ChR2_rep2","FMRP_ChR2_rep3", "FMRP_ChR2_rep4")), na.rm=TRUE),
         Control_max_tags=max(c_across(c("FMRP_Control_rep1", "FMRP_Control_rep2","FMRP_Control_rep3","FMRP_Control_rep4")), na.rm=TRUE))

counts_BC4_tags5=counts_BC4 %>% filter(ChR2_max_tags >= 5 | Control_max_tags >= 5)

nrow(counts_BC4_tags5) #3655

#add in width and gene names
tx_gene_length=mm10_Tx_final_noNAs[,c("transcript_id","gene_id","gene_name","entrez","width")]
counts_BC_genenames=merge(counts_BC4_tags5,tx_gene_length,by.x="Transcript_id", by.y="transcript_id")

counts_BC_genenames=counts_BC_genenames %>% relocate(Transcript_id,
                                                     gene_id,
                                                     gene_name,
                                                     entrez,
                                                     width,
                                                     FMRP_ChR2_rep1,
                                                     FMRP_ChR2_rep2,
                                                     FMRP_ChR2_rep3,
                                                     FMRP_ChR2_rep4,
                                                     FMRP_Control_rep1,
                                                     FMRP_Control_rep2,
                                                     FMRP_Control_rep3,
                                                     FMRP_Control_rep4,
                                                     BC_Control,
                                                     BC_ChR2)
#calculate tags per unit length
counts_forCLIP_TPM <- counts_BC_genenames %>%
  pivot_longer(names_to = "Sample", values_to = "tags", cols = c(6:13)) %>%
  distinct %>%
  mutate(rpk = tags/(width/1000))

#get scaling factors for read normalization
scaling <- counts_forCLIP_TPM %>%
  group_by(Sample) %>%
  summarise(scaling.factor = (sum(rpk))/1000000)

counts_forCLIP_TPM <- left_join(counts_forCLIP_TPM, scaling, by = "Sample")

counts_forCLIP_TPM <- counts_forCLIP_TPM %>%
  mutate(tpm = rpk/scaling.factor) 

#this should equal 1 million for each sample
counts_forCLIP_TPM %>% group_by(Sample) %>% summarise(sum(tpm))

All.FMRP.tags_wider=counts_forCLIP_TPM %>% 
  dplyr::select(Transcript_id,Sample,tags,tpm) %>%
  pivot_wider(names_from = Sample, values_from = c(tags,tpm), values_fill = 0)  

All.FMRP.tags_TPM=merge(All.FMRP.tags_wider, counts_BC_genenames %>% dplyr::select(Transcript_id,BC_ChR2,BC_Control), by.x="Transcript_id", by.y="Transcript_id")

#now merge with RiboTag data
nrow(RiboTag_mm10_gtf_merged_highestExp) #14353
nrow(All.FMRP.tags_TPM) #3655

FMRP_RiboTag_merged <- merge(RiboTag_mm10_gtf_merged_highestExp, All.FMRP.tags_TPM, by.x= "Transcript_ID", by.y= "Transcript_id", all.x=TRUE, all.y=FALSE) 

FMRP_RiboTag_merged <- FMRP_RiboTag_merged %>% 
  mutate_at(c(25:42), ~replace_na(.,0))

FMRP_RiboTag_merged=FMRP_RiboTag_merged %>% 
  ungroup() %>%
  mutate(CLIP_Rest_rep1_TPM=case_when(tpm_FMRP_Control_rep1 > 0 ~ tpm_FMRP_Control_rep1, tpm_FMRP_Control_rep1 <= 0 ~ 1),                                         
         CLIP_Rest_rep2_TPM=case_when(tpm_FMRP_Control_rep2 > 0 ~ tpm_FMRP_Control_rep2, tpm_FMRP_Control_rep2 <= 0 ~ 1),
         CLIP_Rest_rep3_TPM=case_when(tpm_FMRP_Control_rep3 > 0 ~ tpm_FMRP_Control_rep3, tpm_FMRP_Control_rep3 <= 0 ~ 1),
         CLIP_Rest_rep4_TPM=case_when(tpm_FMRP_Control_rep4 > 0 ~ tpm_FMRP_Control_rep4, tpm_FMRP_Control_rep4 <= 0 ~ 1),
         CLIP_Opto_rep1_TPM=case_when(tpm_FMRP_ChR2_rep1 > 0 ~ tpm_FMRP_ChR2_rep1, tpm_FMRP_ChR2_rep1 <= 0 ~ 1),
         CLIP_Opto_rep2_TPM=case_when(tpm_FMRP_ChR2_rep2 > 0 ~ tpm_FMRP_ChR2_rep2, tpm_FMRP_ChR2_rep2 <= 0 ~ 1),
         CLIP_Opto_rep3_TPM=case_when(tpm_FMRP_ChR2_rep3 > 0 ~ tpm_FMRP_ChR2_rep3, tpm_FMRP_ChR2_rep3 <= 0 ~ 1),
         CLIP_Opto_rep4_TPM=case_when(tpm_FMRP_ChR2_rep4 > 0 ~ tpm_FMRP_ChR2_rep4, tpm_FMRP_ChR2_rep4 <= 0 ~ 1))

FMRP_RiboTag_merged= FMRP_RiboTag_merged %>%
  rowwise() %>%
  mutate(avg_rest_CLIP_TPM   = mean(c_across(c("CLIP_Rest_rep1_TPM", "CLIP_Rest_rep2_TPM", "CLIP_Rest_rep3_TPM", "CLIP_Rest_rep4_TPM")), na.rm=TRUE),
         avg_opto_CLIP_TPM   = mean(c_across(c("CLIP_Opto_rep1_TPM", "CLIP_Opto_rep2_TPM", "CLIP_Opto_rep3_TPM", "CLIP_Opto_rep4_TPM")), na.rm=TRUE),
         avg_FC_CLIP_TPM     = avg_opto_CLIP_TPM/avg_rest_CLIP_TPM)


FMRP_RiboTag_merged_log2= FMRP_RiboTag_merged %>% 
  mutate(
    across(
      .cols = c(avg_rest_RT_TPM, avg_opto_RT_TPM, 43:53),
      .fns = log2,
      .names = "{.col}_log2"
    ))

################################################################################################################################################################
#CLIP score generation 
FMRP_RiboTag_merged_both=FMRP_RiboTag_merged_log2 %>% dplyr::filter(avg_rest_RT_TPM_log2 > 0 & avg_opto_RT_TPM_log2 > 0 ) 
nrow(FMRP_RiboTag_merged_both)  #12167

# full Rest clip score code
df <- FMRP_RiboTag_merged_both
my_x <- "avg_rest_RT_TPM_log2"
Treatment <- "Rest"
Rep <- "rep1"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regression line
            add.params = list(color = "red"), # Customize regression line
            cor.coef = TRUE, # Add correlation coefficient
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2A.pdf"), plot = myplot_pearson)

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
myplot_loess
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2B.pdf"), plot = myplot_pearson)
fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit

myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Rep <- "rep3"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2C.pdf"), plot = myplot_pearson)
fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2D.pdf"), plot = myplot_pearson)
fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

################################################################################
# full Opto clip score code
my_x <- "avg_opto_RT_TPM_log2"
Treatment <- "Opto"
Rep <- "rep1"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2E.pdf"), plot = myplot_pearson)

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep2
Rep <- "rep2"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2F.pdf"), plot = myplot_pearson)

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep3
Rep <- "rep3"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2G.pdf"), plot = myplot_pearson)

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

#rep4
Rep <- "rep4"
my_y <- paste0("CLIP_", Treatment, "_", Rep, "_TPM_log2")
my_y

myplot_pearson <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggscatter(., x= my_x, y= my_y,
            color = "black", shape = 19, size = .3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "red"), # Customize reg. line
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
            xlab=(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")),
            ylab=(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM")))
myplot_pearson
ggsave(paste0(Outdirectory, "/CLIP_", Treatment, "_", Rep, "_CLIPscores_pearson_S2H.pdf"), plot = myplot_pearson)

fm <- as.formula(paste(my_y, "~", my_x))
fit <- lm(fm, data = df, subset = (get(my_y) >= 0))
fit
myplot_loess <- df %>%
  dplyr::filter(get(my_y) >= 0) %>%
  ggplot() +
  aes_string(
    x = my_x,
    y = my_y
  ) +
  geom_point(cex = 0.3, pch = 19) +
  geom_abline(
    slope = coef(fit)[[2]],
    intercept = coef(fit)[["(Intercept)"]], color = "red"
  ) +
  xlab(paste0("log2 ", Treatment, " Camk2a RiboTag TPM")) +
  ylab(paste0("log2 ", Treatment, "_", Rep, " Camk2a Fmr1cTag CLIP TPM"))
newcolumn <- paste0(Treatment, "_", Rep, "_CLIPscore")
newcolumn
df <- df %>% mutate(!!newcolumn := (get(my_y) - (coef(fit)[[2]] * get(my_x)) + coef(fit)[["(Intercept)"]]))
colnames(df)

CLIPscores_final=df

################################################################################################################################################################
CLIPscores_merged=CLIPscores_final %>% 
  rowwise() %>%
  mutate(avg_rest_CLIPscores = mean(c_across(c("Rest_rep1_CLIPscore","Rest_rep2_CLIPscore","Rest_rep3_CLIPscore","Rest_rep4_CLIPscore")), na.rm=TRUE),
         avg_opto_CLIPscores = mean(c_across(c("Opto_rep1_CLIPscore","Opto_rep2_CLIPscore","Opto_rep3_CLIPscore","Opto_rep4_CLIPscore")), na.rm=TRUE))

CLIPscores_merged <- CLIPscores_merged %>% 
  rowwise() %>%
  mutate(color_rest_CLIPtpm = case_when(avg_rest_CLIPscores > 2 ~ "Stringent_Targets_Rest",
                                        avg_rest_CLIPscores <= 2 & avg_rest_CLIPscores > 1 ~ "High_binding_Rest",
                                        avg_rest_CLIPscores <=1 & avg_rest_CLIPscores > 0 ~ "Low_binding_Rest",
                                        avg_rest_CLIPscores <=0 ~ "NonFMRP_target")) %>%
  mutate(color_opto_CLIPtpm = case_when(avg_opto_CLIPscores > 2 ~ "Stringent_Targets_Opto",
                                        avg_opto_CLIPscores <= 2 & avg_opto_CLIPscores > 1 ~ "High_binding_Opto",
                                        avg_opto_CLIPscores <=1 & avg_opto_CLIPscores > 0 ~ "Low_binding_Opto",
                                        avg_opto_CLIPscores <=0 ~ "NonFMRP_target"))

nrow(subset(CLIPscores_merged, color_rest_CLIPtpm == "Stringent_Targets_Rest")) #1602
nrow(subset(CLIPscores_merged, color_rest_CLIPtpm == "High_binding_Rest"))  #915
nrow(subset(CLIPscores_merged, color_rest_CLIPtpm == "Low_binding_Rest"))   #580
nrow(subset(CLIPscores_merged, color_rest_CLIPtpm == "NonFMRP_target"))     #9070
nrow(subset(CLIPscores_merged, avg_rest_CLIP_TPM_log2 <= 0))                #8569

nrow(subset(CLIPscores_merged, color_opto_CLIPtpm == "Stringent_Targets_Opto")) #1027 
nrow(subset(CLIPscores_merged, color_opto_CLIPtpm == "High_binding_Opto"))      #596
nrow(subset(CLIPscores_merged, color_opto_CLIPtpm == "Low_binding_Opto"))       #717
nrow(subset(CLIPscores_merged, color_opto_CLIPtpm == "NonFMRP_target"))         #9827
nrow(subset(CLIPscores_merged, avg_opto_CLIP_TPM_log2 <= 0))                    #8596

Rest_avg_CLIPscore_plot=ggplot(CLIPscores_merged %>% filter(color_rest_CLIPtpm == "NonFMRP_target" )) +
  geom_point(aes(x=avg_rest_RT_TPM_log2, y=log2(avg_rest_CLIP_TPM), color="NonFMRP_target"), size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(color_rest_CLIPtpm == "Low_binding_Rest"), aes(x=avg_rest_RT_TPM_log2, y=log2(avg_rest_CLIP_TPM), color="Low_binding_Rest"), size =2 ) +
  geom_point(data= CLIPscores_merged %>% filter(color_rest_CLIPtpm == "High_binding_Rest"), aes(x=avg_rest_RT_TPM_log2, y=log2(avg_rest_CLIP_TPM), color="High_binding_Rest"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(color_rest_CLIPtpm == "Stringent_Targets_Rest"), aes(x=avg_rest_RT_TPM_log2, y=log2(avg_rest_CLIP_TPM), color="Stringent_Targets_Rest"), size =2) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Low_binding_Rest = "grey",
                                High_binding_Rest = "lightskyblue1",
                                Stringent_Targets_Rest = "dodgerblue")) + 
  xlab("Rest RiboTag TPM (log2)") +
  ylab("Rest FMRP CLIP TPM (log2FC)") +
  xlim(c(0,15)) +
  ylim(c(0,17)) +
  theme_classic() 
Rest_avg_CLIPscore_plot
ggsave(filename =paste0(Outdirectory,"/Figure6A_Control_avg_CLIPscore_plot.pdf"), plot =Rest_avg_CLIPscore_plot)

Opto_avg_CLIPscore_plot=ggplot(CLIPscores_merged %>% filter(color_opto_CLIPtpm == "NonFMRP_target" )) +
  geom_point(aes(x=avg_opto_RT_TPM_log2, y=log2(avg_opto_CLIP_TPM), color="NonFMRP_target"),size =2) + 
  geom_point(data= CLIPscores_merged %>% filter(color_opto_CLIPtpm == "Low_binding_Opto"), aes(x=avg_opto_RT_TPM_log2, y=log2(avg_opto_CLIP_TPM), color="Low_binding_Opto"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(color_opto_CLIPtpm == "High_binding_Opto"), aes(x=avg_opto_RT_TPM_log2, y=log2(avg_opto_CLIP_TPM), color="High_binding_Opto"), size =2) +
  geom_point(data= CLIPscores_merged %>% filter(color_opto_CLIPtpm == "Stringent_Targets_Opto"), aes(x=avg_opto_RT_TPM_log2, y=log2(avg_opto_CLIP_TPM), color="Stringent_Targets_Opto"), size =2) +
  scale_color_manual(name = "FMRP CLIP score", 
                     values = c(NonFMRP_target ="black",
                                Low_binding_Opto = "grey",
                                High_binding_Opto = "lightcoral",
                                Stringent_Targets_Opto = "red")) + 
  xlab("Opto RiboTag TPM (log2)") +
  ylab("Opto FMRP CLIP TPM (log2FC)") +
  xlim(c(0,15)) +
  ylim(c(0,17)) +
  theme_classic()
Opto_avg_CLIPscore_plot
ggsave(filename =paste0(Outdirectory,"/Figure6B_Opto_avg_CLIPscore_plot.pdf"), plot =Opto_avg_CLIPscore_plot)

CLIPscores_merged <- CLIPscores_merged %>%
  rowwise() %>%
  mutate(color_CLIPavg = case_when(avg_rest_CLIPscores > 2 & avg_opto_CLIPscores > 2 ~ "Stringent_Targets_Both",
                                   avg_rest_CLIPscores > 2 & avg_opto_CLIPscores <= 2 ~ "Stringent_Targets_Rest",
                                   avg_rest_CLIPscores <= 2 & avg_opto_CLIPscores > 2 ~ "Stringent_Targets_Opto",
                                   avg_rest_CLIPscores <= 2 & avg_opto_CLIPscores <= 2 ~ "NonFMRP_target"))

restvsopto_avg=ggplot(CLIPscores_merged, aes(x=avg_rest_CLIPscores, y=avg_opto_CLIPscores, color=color_CLIPavg)) + 
  scale_color_manual(name = "FMRP CLIP score", values = c(Stringent_Targets_Both = "green4", Stringent_Targets_Rest = "dodgerblue", Stringent_Targets_Opto = "red3", NonFMRP_target = "darkgray")) +
  geom_point() +
  xlim(c(-8,15)) +
  ylim(c(-8,15)) +
  theme_classic() +
  geom_hline(yintercept = 2, colour = "grey", linewidth=.3) + 
  geom_vline(xintercept = 2, colour = "grey", linewidth=.3) 
restvsopto_avg
ggsave(filename =paste0(Outdirectory,"/Figure6C_Control_vs_Opto_avgCLIPscores.pdf"), plot =restvsopto_avg)

nrow(subset(CLIPscores_merged, color_CLIPavg == "Stringent_Targets_Opto")) #129
nrow(subset(CLIPscores_merged, color_CLIPavg == "Stringent_Targets_Rest")) #704
nrow(subset(CLIPscores_merged, color_CLIPavg == "Stringent_Targets_Both")) #898
nrow(subset(CLIPscores_merged, color_CLIPavg == "NonFMRP_target")) #10436
 
################################################################################################################################################################
CLIPscores_matrix_raw=CLIPscores_merged %>% dplyr::select(gene_name,
                                                          Rest_rep1_CLIPscore,
                                                          Rest_rep2_CLIPscore,
                                                          Rest_rep3_CLIPscore,
                                                          Rest_rep4_CLIPscore,
                                                          Opto_rep1_CLIPscore,
                                                          Opto_rep2_CLIPscore,
                                                          Opto_rep3_CLIPscore,
                                                          Opto_rep4_CLIPscore)

nrow(CLIPscores_matrix_raw) #12167

CLIPscores_matrix_zeros <- CLIPscores_matrix_raw %>% 
  mutate_at(c(2:9), ~ifelse(. < 0, 0, .))

CLIPscores_matrix_use=CLIPscores_matrix_zeros %>% column_to_rownames(var="gene_name")

Limma_CLIPscores <- DGEList(CLIPscores_matrix_use)
Limma_CLIPscores <- calcNormFactors(Limma_CLIPscores)

group=as.factor(c("Rest", "Rest", "Rest","Rest",
                  "Opto","Opto","Opto","Opto"))
Limma_CLIPscores$samples$group <- group

Limma_CLIPscores

mm <- model.matrix(~0 + group)
y <- voom(Limma_CLIPscores, mm, plot = T)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupOpto - groupRest, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
colnames(top.table)=paste(colnames(top.table), "limma", sep = ".")
top.table <- tibble::rownames_to_column(top.table, "gene_name")

CLIPscores_Limma=merge(CLIPscores_merged,top.table,by.x="gene_name",by.y="gene_name")

CLIPscores_Limma=CLIPscores_Limma[order(CLIPscores_Limma$P.Value.limma), ]

vol_plot=ggplot(CLIPscores_Limma) +
  geom_point(data=CLIPscores_Limma %>% filter(P.Value.limma >= 0.05),aes(x = logFC.limma, y = -log10(P.Value.limma)), color="grey") +
  geom_point(data=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma > 0),aes(x = logFC.limma, y = -log10(P.Value.limma)), color="red3") +
  geom_point(data=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma < 0),aes(x = logFC.limma, y = -log10(P.Value.limma)), color="dodgerblue") +
  xlab(paste("Log2 Fold Change\n","Opto vs Rest CLIP scores", sep=" ")) + 
  ylab(expression(-log[10]("P-value"))) + 
  geom_hline(yintercept = 1.3, linetype='dotted', colour = "darkgrey") +
  theme_classic() +
  scale_y_continuous(trans = "log1p") + 
  xlim(-3.5,3.5) +
  theme(text=element_text(size=16), legend.position = 'none') +
  annotate("text", 
           x = .5, 
           y = 7.5, 
           label = paste0("Decreased (", nrow(CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma < 0)),")"),
           color = "dodgerblue") +
  annotate("text", 
           x = .5, 
           y = 8.5, 
           label = paste0("Increased (", nrow(CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma > 0)),")"),
           color = "red3") 
vol_plot
ggsave(filename = paste0(Outdirectory,"/","Figure6D_FMRP_CLIPscores_limma_Volcanoplot_unlabeled.pdf"), plot =vol_plot)

vol_plot_labeled_text=vol_plot + 
  geom_text_repel(data=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma > 0) %>% arrange(P.Value.limma) %>% dplyr::slice_head(n=25), aes(x=logFC.limma, y=-log10(P.Value.limma),label=gene_name),max.overlaps=1000, show.legend=FALSE) +
  geom_text_repel(data=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & logFC.limma < 0) %>% arrange(P.Value.limma) %>% dplyr::slice_head(n=25), aes(x=logFC.limma, y=-log10(P.Value.limma),label=gene_name),max.overlaps=1000, show.legend=FALSE) 
vol_plot_labeled_text
ggsave(filename = paste0(Outdirectory,"/","Figure6D_FMRP_CLIPscores_limma_Volcanoplot_labeled_text.pdf"), plot =vol_plot_labeled_text)

###GSEA
CLIPscores_Limma_sig=CLIPscores_Limma %>% dplyr::filter(P.Value.limma < 0.05)
original_gene_list <- CLIPscores_Limma_sig$t.limma

# name the vector
names(original_gene_list) <- CLIPscores_Limma_sig$gene_name

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = Inf, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")

gse_result=gse@result

gsea_res.sig <- gse_result %>%
  filter(pvalue<0.05)

gsea_res.sig.opto <- gsea_res.sig %>%
  filter(NES > 0)

gsea_res.sig.opto %>% arrange(desc(NES)) %>% slice(1:40)

gsea_res.sig.rest <- gsea_res.sig %>%
  filter(NES < 0)

gsea_res.sig.rest %>% arrange(NES) %>% slice(1:40)

categories=c("immune response",
             "transcription factor binding",
             "negative regulation of inflammatory response",
             "growth factor receptor binding",
             "chromatin binding",
             "RNA localization",
             "transcription by RNA polymerase II",
             "regulation of peptidyl-serine phosphorylation",
             "regulation of AMPA receptor activity",
             "mitochondrial fusion",
             "postsynaptic density membrane",
             "glutamatergic synapse",
             "intracellular calcium ion homeostasis",
             "Schaffer collateral - CA1 synapse")

categories_orderedbyNES=gsea_res.sig %>% filter(Description %in% categories) %>% arrange(NES) %>% select(Description)

clipscore_limma_gsea=ggplot(subset(gsea_res.sig, Description %in% categories), aes(NES, Description, color=p.adjust, size= setSize)) +
  geom_point() +
  scale_color_viridis_c() +
  scale_size_binned()+
  scale_y_discrete(limits=factor(categories_orderedbyNES$Description)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) 
clipscore_limma_gsea
ggsave(filename = paste0(Outdirectory,"/Figure6E_Clipscore_Limma_GSEA.pdf"), plot =clipscore_limma_gsea, width = 7, height=8, unit="in")

##############################################################################################################################
Cellbody_v_Neurites_CLIPscores <- read.csv(file.path(Datadirectory,"Hale_etal_2019_PMID_34939924_FMRP_CBvD_CLIPscores.csv"))

OptoCLIP_CBvN_CLIPscores_merged=merge(CLIPscores_Limma,
                                      Cellbody_v_Neurites_CLIPscores, by.x="gene_name", by.y="GeneName")

OptoCLIP_CBvN_CLIPscores_merged <- OptoCLIP_CBvN_CLIPscores_merged %>%
  rowwise() %>%
  mutate(color_CLIPscores_Limma = case_when(P.Value.limma >= 0.05 ~ "NotDiffBound",
                                         P.Value.limma < 0.05 & logFC.limma < 0 ~ "CLIP_Rest",
                                         P.Value.limma < 0.05 & logFC.limma > 0 ~ "CLIP_Opto")) 

#the following violin plots show that stats comparing based the comparisons I set below 

my_comparisons <- list( c("NotDiffBound", "CLIP_Rest"), c("NotDiffBound", "CLIP_Opto"), c("CLIP_Rest", "CLIP_Opto") )

violin_CB=ggviolin(OptoCLIP_CBvN_CLIPscores_merged, x = "color_CLIPscores_Limma", y = "CB.clip.score", fill = "color_CLIPscores_Limma",
                   palette = c("grey", "dodgerblue", "red3"),
                   add = "boxplot", add.params = list(fill = "white"),
                   ylab="FMRP CLIP Score Cell Bodies",
                   xlab="")+
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test") +
  theme_classic() +
  ylim(-5,10) +
  theme(text=element_text(size=16), legend.position = 'none') 
violin_CB
ggsave(filename = paste0(Outdirectory,"/Figure6F_OptoCLIP_CB_Violinplot.pdf"), plot =violin_CB, width = 6, height=6, unit="in")

violin_Neurites=ggviolin(OptoCLIP_CBvN_CLIPscores_merged, x = "color_CLIPscores_Limma", y = "dendritic.CLIP.score", fill = "color_CLIPscores_Limma",
                         palette = c("grey", "dodgerblue", "red3"),
                         add = "boxplot", add.params = list(fill = "white"),
                         ylab="FMRP CLIP Score Neurites",
                         xlab="")+
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test") +
  theme_classic() +
  ylim(-5,10) +
  theme(text=element_text(size=16), legend.position = 'none') 
violin_Neurites
ggsave(filename = paste0(Outdirectory,"/Figure6G_OptoCLIP_Neurites_Violinplot.pdf"), plot =violin_Neurites, width = 6, height=6, unit="in")

################################################################################################################################
#combining RiboTag data with CLIP scores

CLIPscores_Limma= CLIPscores_Limma %>% ungroup()

CLIP_RiboTag_plot=ggplot(CLIPscores_Limma) +   
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue > 0.05 & logFC.limma < 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_NS"), size =2 ) + 
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma < 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_Opto"), size =2 ) + 
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma < 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_Rest"), size =2 ) + 
  ylab("RiboTag TPM Opto vs Rest log2FC") +
  xlab("FMRP CLIP Scores Opto vs Rest log2FC") +
  theme_classic() +
  xlim(c(-3.5,0)) + 
  ylim(c(-3.5,3.5)) +
  geom_hline(yintercept = 0, colour = "black", linewidth=.3) + 
  scale_color_manual(name = "RiboTag & FMRP CLIP Score", 
                     values = c(CLIP_Rest_Ribo_NS = "dodgerblue",
                                CLIP_Rest_Ribo_Opto = "violet",
                                CLIP_Rest_Ribo_Rest = "darkorange")) + 
  geom_label_repel(data=CLIPscores_Limma %>% dplyr::filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma < 0) %>% arrange(desc(log2FoldChange)) %>% dplyr::slice(1:25), aes(x=logFC.limma, y=log2FoldChange, label=gene_name),nudge_x = -.5,nudge_y=.5,color="violet", segment.color="black",max.overlaps=Inf, show.legend=FALSE) +
  geom_label_repel(data=CLIPscores_Limma %>% dplyr::filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma < 0) %>% arrange(log2FoldChange) %>% dplyr::slice(1:25), aes(x=logFC.limma, y=log2FoldChange, label=gene_name),nudge_x = -.5,nudge_y=-.5,color="darkorange", segment.color="black",max.overlaps=Inf, show.legend=FALSE) 
CLIP_RiboTag_plot
ggsave(filename = paste0(Outdirectory,"/","Figure7A_CLIP_RiboTag_plot.pdf"), plot =CLIP_RiboTag_plot)

#GO analysis
CLIP_Rest_Ribo_Opto=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma < 0)
CLIP_Rest_Ribo_Rest=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma < 0)

nrow(CLIP_Rest_Ribo_Opto) #94
nrow(CLIP_Rest_Ribo_Rest) #72

CLIP_Rest_Ribo_Opto$gene_name #"Gls" "Fndc4"         "Zfp512b"       "Cnr1"          "Pou6f1"        "Strbp"         "Sh3bp4"        "2510009E07Rik" "Usp34"         "Pcdh10"     "Dyrk2"         "Wdr6"          "Usp5"          "Kifc2"         "Plxnd1"        "Pou3f3"        "Bltp2"         "Hunk"             
CLIP_Rest_Ribo_Rest$gene_name #"Csrnp1" "Nfe2l2"   "Zc3hav1"  "Aldoc"    "Ahnak"    "Tjp1"     "Zfp36l1"  "Glul"     "Zfp36"    "Sik1"     "Ccn1"     "Clu"      "Ifrd1"    "Atp1a2"   "Ppp1r15b" "Ddx17"    "Ifit3b"   "Bahcc1"   "Ptgs2"   

GO_CLIP_Rest_Ribo_Opto= enrichGO(gene = CLIP_Rest_Ribo_Opto$entrez,
                        OrgDb = 'org.Mm.eg.db', 
                        keyType = "ENTREZID",
                        readable = T,
                        ont = "all",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
GO_CLIP_Rest_Ribo_Opto_results=GO_CLIP_Rest_Ribo_Opto@result

CLIP_Rest_Ribo_Opto_terms_ordered=c("axon guidance",
                                    "synaptic transmission, glutamatergic",
                                    "microtubule binding",
                                    "cysteine-type peptidase activity",
                                    "synaptic vesicle localization",
                                    "miRNA binding")

CLIP_Rest_Ribo_Opto_plot=ggplot(GO_CLIP_Rest_Ribo_Opto_results %>% dplyr::filter(Description %in% CLIP_Rest_Ribo_Opto_terms_ordered), aes(Count, Description)) +
  geom_bar(stat = "identity", fill="violet") +
  scale_y_discrete(limits=factor(CLIP_Rest_Ribo_Opto_terms_ordered)) +
  theme_classic()
CLIP_Rest_Ribo_Opto_plot
ggsave(filename =paste0(Outdirectory,"/Figure7B_CLIP_Rest_Ribo_Opto_GOterms.pdf"), width = 8, height = 8, units = "in", plot =CLIP_Rest_Ribo_Opto_plot)

GO_CLIP_Rest_Ribo_Rest=enrichGO(gene = CLIP_Rest_Ribo_Rest$entrez,
                                 OrgDb = 'org.Mm.eg.db', 
                                 keyType = "ENTREZID",
                                 readable = T,
                                 ont = "all",
                                 pvalueCutoff = 0.05, 
                                 qvalueCutoff = 0.10)
GO_CLIP_Rest_Ribo_Rest_results=GO_CLIP_Rest_Ribo_Rest@result

CLIP_Rest_Ribo_Rest_terms_ordered=c("cytoplasmic ribonucleoprotein granule",
                                    "regulation of translation",
                                    "mRNA binding",
                                    "positive regulation of mRNA catabolic process",
                                    "14-3-3 protein binding",
                                    "3'-UTR-mediated mRNA destabilization")

CLIP_Rest_Ribo_Rest_plot=ggplot(GO_CLIP_Rest_Ribo_Rest_results %>% dplyr::filter(Description %in% CLIP_Rest_Ribo_Rest_terms_ordered), aes(Count, Description)) +
  geom_bar(stat = "identity", fill="darkorange") +
  scale_y_discrete(limits=factor(CLIP_Rest_Ribo_Rest_terms_ordered)) +
  theme_classic()
CLIP_Rest_Ribo_Rest_plot
ggsave(filename =paste0(Outdirectory,"/Figure7C_CLIP_Rest_Ribo_Rest_GOterms.pdf"), width = 8, height = 8, units = "in", plot =CLIP_Rest_Ribo_Rest_plot)

################################################################################################################################
#combining RiboTag data with CLIP scores for CLIP Up Transcripts
#FigureS3A

CLIPscores_Limma= CLIPscores_Limma %>% ungroup()

CLIP_RiboTag_plot_CLIP_Up=ggplot(CLIPscores_Limma) +   
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue > 0.05 & logFC.limma > 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_NS"), size =2 ) + 
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma > 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_Opto"), size =2 ) + 
  geom_point(data= CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma > 0), aes(x=logFC.limma, y=log2FoldChange, color="CLIP_Rest_Ribo_Rest"), size =2 ) + 
  ylab("RiboTag TPM Opto vs Rest log2FC") +
  xlab("FMRP CLIP Scores Opto vs Rest log2FC") +
  theme_classic() +
  ylim(c(-1,1)) +
  xlim(c(0,2)) +
  geom_hline(yintercept = 0, colour = "black", linewidth=.3) + 
  scale_color_manual(name = "RiboTag & FMRP CLIP Score", 
                     values = c(CLIP_Rest_Ribo_NS = "red3",
                                CLIP_Rest_Ribo_Opto = "green3",
                                CLIP_Rest_Ribo_Rest = "darkorchid2")) +
  geom_label_repel(data=CLIPscores_Limma %>% dplyr::filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma > 0) %>% arrange(desc(log2FoldChange)) %>% dplyr::slice(1:25), aes(x=logFC.limma, y=log2FoldChange, label=gene_name),nudge_x = -.1,nudge_y=.1,color="green3", segment.color="black",max.overlaps=Inf, show.legend=FALSE) +
  geom_label_repel(data=CLIPscores_Limma %>% dplyr::filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma > 0) %>% arrange(log2FoldChange) %>% dplyr::slice(1:25), aes(x=logFC.limma, y=log2FoldChange, label=gene_name),nudge_x = -.1,nudge_y=-.1,color="darkorchid2", segment.color="black",max.overlaps=Inf, show.legend=FALSE) 
CLIP_RiboTag_plot_CLIP_Up
ggsave(filename = paste0(Outdirectory,"/","FigureS3A_CLIP_Up_RiboTag_plot.pdf"), plot =CLIP_RiboTag_plot_CLIP_Up)

CLIP_Opto_Ribo_Opto=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange > 0 & logFC.limma > 0)
CLIP_Opto_Ribo_Rest=CLIPscores_Limma %>% filter(P.Value.limma < 0.05 & pvalue < 0.05 & log2FoldChange < 0 & logFC.limma > 0)

nrow(CLIP_Opto_Ribo_Opto) #9
nrow(CLIP_Opto_Ribo_Rest) #10

CLIP_Opto_Ribo_Opto$gene_name #"Gnai1"    "Prkaca"   "Syt4"     "Mir9-2hg" "Ctxnd1"   "Fmn1"     "Gnaq"     "Mink1"    "Cenpb"             
CLIP_Opto_Ribo_Rest$gene_name #"Synpo"   "Lnpk"    "Ahcyl2"  "Akirin1" "Srp54c"  "Dnajc5"  "Slc2a3"  "Acvr2a"  "Zzz3"    "Stat3"  
