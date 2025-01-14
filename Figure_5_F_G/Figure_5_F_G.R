library(rtracklayer)
library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tximport)
library(DESeq2)
library(clusterProfiler)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","ruthasinger_github","Opto-CLIP_v2","Figure_5_F_G")
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

RiboTag_Metadata <- mutate_at(RiboTag_Metadata, vars(Group), as.factor)

RT_dds <- DESeqDataSetFromTximport(RiboTag_txi.tx, RiboTag_Metadata, ~ Group)

keep <- rowSums(counts(RT_dds)) >= 1
RT_dds <- RT_dds[keep, ]
RT_dds <- DESeq(RT_dds)

res=results(RT_dds, contrast=c("Group", "ChR2","Control"))
res=res[order(res$padj), ]
resdata=as.data.frame(res)
resdata <- tibble::rownames_to_column(resdata, "Transcript_ID")

RiboTag_salmon_TPM=as.data.frame(RiboTag_txi.tx$abundance)
RiboTag_salmon_TPM <- tibble::rownames_to_column(RiboTag_salmon_TPM, "trancript_ID")
colnames(RiboTag_salmon_TPM) <- c("transcript_ID",paste(colnames(RiboTag_salmon_TPM[2:ncol(RiboTag_salmon_TPM)]), "TPM", sep = "."))

RiboTag_resdata_TPM=merge(resdata,RiboTag_salmon_TPM, by.x="Transcript_ID", by.y="transcript_ID", all.x=TRUE, all.y=FALSE)

#add in gtf info
RiboTag_mm10_gtf_merged <- merge(RiboTag_resdata_TPM, mm10_Tx_final_noNAs, by.x="Transcript_ID", by.y="transcript_id")

#should be true
length(unique(RiboTag_mm10_gtf_merged$Transcript_ID)) == nrow(RiboTag_mm10_gtf_merged)

#should be false
length(unique(RiboTag_mm10_gtf_merged$gene_name)) == nrow(RiboTag_mm10_gtf_merged)

nrow(RiboTag_mm10_gtf_merged) #60270

colnames(RiboTag_mm10_gtf_merged)

RiboTag_mm10_gtf_merged_exp=RiboTag_mm10_gtf_merged %>% 
  dplyr::filter(RiboTag_ChR2_30min_IP_rep1_salmon.TPM > 0 & 
                RiboTag_ChR2_30min_IP_rep2_salmon.TPM > 0 & 
                RiboTag_ChR2_30min_IP_rep3_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep1_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep2_salmon.TPM > 0 & 
                RiboTag_Control_30min_IP_rep3_salmon.TPM > 0)

nrow(RiboTag_mm10_gtf_merged_exp) #30157

#Figure 5F
RiboTag_mm10_gtf_merged_exp=RiboTag_mm10_gtf_merged_exp[order(RiboTag_mm10_gtf_merged_exp$padj), ]

vol_plot=ggplot(RiboTag_mm10_gtf_merged_exp) +
  geom_point(data=RiboTag_mm10_gtf_merged_exp %>% filter(padj >= 0.05),aes(x = log2FoldChange, y = -log10(padj)), color="grey") +
  geom_point(data=RiboTag_mm10_gtf_merged_exp %>% filter(padj < 0.05 & log2FoldChange > 0),aes(x = log2FoldChange, y = -log10(padj)), color="red3") +
  geom_point(data=RiboTag_mm10_gtf_merged_exp %>% filter(padj < 0.05 & log2FoldChange < 0),aes(x = log2FoldChange, y = -log10(padj)), color="dodgerblue") +
  xlab(paste("Log2 Fold Change\n","Opto vs Control RiboTag", sep=" ")) + 
  ylab(expression(-log[10]("P-adjusted"))) + 
  geom_hline(yintercept = 1.3, linetype='dotted', colour = "darkgrey") +
  theme_classic() +
  scale_y_continuous(trans = "log1p") + 
  xlim(-5,5) +
  theme(text=element_text(size=16), legend.position = 'none') +
  annotate("text", 
           x = 3, 
           y = 30, 
           label = paste0("Decreased (", nrow(RiboTag_mm10_gtf_merged_exp %>% filter(padj < 0.05 & log2FoldChange < 0)),")"),
           color = "dodgerblue") +
  annotate("text", 
           x = 3, 
           y = 25, 
           label = paste0("Increased (", nrow(RiboTag_mm10_gtf_merged_exp %>% filter(padj < 0.05 & log2FoldChange > 0)),")"),
           color = "red3") 
vol_plot
ggsave(filename = paste0(Outdirectory,"/Figure5F.pdf"), plot =vol_plot)

################################################################################################################################
#Figure 5G- GSEA
RT_forGO_Up=subset(RiboTag_mm10_gtf_merged_exp, padj < 0.05 & log2FoldChange > 0 & is.na(entrez) == FALSE)
nrow(RT_forGO_Up) #365

df=RT_forGO_Up
go_enrich_all= enrichGO(gene = df$entrez,
                        OrgDb = 'org.Mm.eg.db', 
                        keyType = "ENTREZID",
                        readable = T,
                        ont = "all",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
godata_RT_forGO_Up=go_enrich_all@result
rownames(godata_RT_forGO_Up)=NULL

godata_RT_forGO_Up %>% arrange(desc(Count)) %>% slice(1:100) %>% select(Description)

RT_terms_Up=c("postsynaptic specialization",
              "neuron to neuron synapse",
              "protein serine/threonine kinase activity",
              "postsynaptic density",
              "GTPase regulator activity",
              "microtubule",
              "cognition",
              "learning or memory")

godata_RT_forGO_Up_use=subset(godata_RT_forGO_Up, Description %in% RT_terms_Up)

RT_forGO_Down=subset(RiboTag_mm10_gtf_merged_exp, padj < 0.05 & log2FoldChange < 0 & is.na(entrez) == FALSE)
nrow(RT_forGO_Down) #418

df=RT_forGO_Down
go_enrich_all= enrichGO(gene = df$entrez,
                        OrgDb = 'org.Mm.eg.db', 
                        keyType = "ENTREZID",
                        readable = T,
                        ont = "all",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
godata_RT_forGO_Down=go_enrich_all@result
rownames(godata_RT_forGO_Down)=NULL

godata_RT_forGO_Down %>% arrange(desc(Count)) %>% slice(1:100) %>% select(Description)

RT_terms_Down=c("regulation of innate immune response",
                "negative regulation of response to external stimulus",
                "regulation of neurogenesis",
                "ubiquitin-like protein ligase binding",
                "protein localization to nucleus",
                "mRNA binding",
                "mRNA processing",
                "RNA splicing")

godata_RT_forGO_Down_use=subset(godata_RT_forGO_Down, Description %in% RT_terms_Down)

RT_GOterms_list=list(godata_RT_forGO_Up_use,godata_RT_forGO_Down_use)
names(RT_GOterms_list)=c("RT_forGO_Up","RT_forGO_Down")

Allgoterms_RiboTag=do.call("rbind", RT_GOterms_list)

Allgoterms_RiboTag$sourcefile=rownames(Allgoterms_RiboTag)
row.names(Allgoterms_RiboTag) <- NULL

Allgoterms_RiboTag= Allgoterms_RiboTag %>% separate(sourcefile, into = c("Sample","two"), sep = "\\.", remove = TRUE) 
Allgoterms_RiboTag= Allgoterms_RiboTag %>% dplyr::select(!two)

RT_GO_plot=ggplot(Allgoterms_RiboTag, aes(Count, Description, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits=c("postsynaptic specialization",
                            "neuron to neuron synapse",
                            "protein serine/threonine kinase activity",
                            "postsynaptic density",
                            "GTPase regulator activity",
                            "microtubule",
                            "cognition",
                            "learning or memory",
                            "regulation of innate immune response",
                            "negative regulation of response to external stimulus",
                            "regulation of neurogenesis",
                            "ubiquitin-like protein ligase binding",
                            "protein localization to nucleus",
                            "mRNA binding",
                            "mRNA processing",
                            "RNA splicing")) +
  theme_classic() +
  scale_fill_manual(name = "RiboTag", 
                    values = c(RT_forGO_Up ="red3",
                               RT_forGO_Down = "dodgerblue"))
RT_GO_plot
ggsave(paste0(Outdirectory, "/Figure5G.pdf"), plot = RT_GO_plot, device = "pdf")
