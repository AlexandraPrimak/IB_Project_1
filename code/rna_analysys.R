library(readxl)
library(data.table)
library(ggplot2)
library(dplyr)
library(ontologyIndex)
library(ontologySimilarity)
library(ggthemes)

# load HPO ontology
data(hpo)


# read gene to HPO mapping: contains frequency of HPO per gene
gene_hpo <- fread("/Users/alexandraprimak/Desktop/genes_to_phenotype.txt")
str(gene_hpo)

hpo_gene <- fread("/Users/alexandraprimak/Desktop/phenotype_to_genes.txt")
colnames(hpo_gene) <- c("HPO_ID", "HPO_name", "entrezID", "geneID", "Additional_Info",  "source", "disease-ID")
hpo_gene <- hpo_gene[ , c("HPO_ID", "HPO_name",  "geneID")]
hpo_gene <- hpo_gene[!duplicated(hpo_gene)]



# read GTEX proteomics data for RNA
annotation_gtex <- read_excel("/Users/alexandraprimak/Desktop/Table_S2_gene_info_at_RNA_level.xlsx" ,
                              sheet = 1 ) %>% as.data.table()

# TS scores
ts_scores <- read_excel("/Users/alexandraprimak/Desktop/Table_S2_gene_info_at_RNA_level.xlsx" ,  
                        sheet = 5 ) %>% as.data.table
ts_scores <- ts_scores[ , 4 : ncol(ts_scores)] #subset necessary columns 
str(ts_scores)


tss <- melt(ts_scores, id.vars = "hgnc_symbol") # reshape data 
setnames(tss, c("geneID", "tissue", "TSscore"))
str(tss) # TS score is a character - wrong

tss[, c("TSscore_all", "NA_value") := tstrsplit(TSscore, ";", fixed=TRUE) ] # split TS score by ; symbol

tss$TSscore_all <- as.numeric(tss$TSscore_all) # convert to numeric
tss$TSscore <- tss$TSscore_all #reassign
str(tss)
tss[ NA_value == "NA_one_rep_in_raw", TSscore := NA] # TSscore has missing value now if low confidence, TSscore_all - contains all of the data


# Create annotation of scores
# TS score reaches 2.5 (SD from the mean of the population distribution) in at least one tissue
# Similarly, if the TS score of a protein is greater than 4 in a tissue and is at least 1.5 higher than the protein’s TS scores in any other tissue, this protein is considered tissue specific.

# Read already prepared annotation
enrichment_results <- read_excel("/Users/alexandraprimak/Desktop/Table_S3_protein_RNA_enrichment_and_concordance_comparison.xlsx" ,  
                                 sheet = 2 ) %>% as.data.table
enrichment_results <- enrichment_results[, c("hgnc_symbol...4", "prt_ench_category", "rna_ench_category" )]
setnames(enrichment_results, c("hgnc_symbol...4", "prt_ench_category", "rna_ench_category" ),
         c("geneID", "Protein_category", "RNA_category"))

enrichment_results_RNA <- enrichment_results[geneID %in% unique(tss$geneID) ,  c("geneID", "RNA_category") ]
enrichment_results_RNA <- enrichment_results_RNA[!duplicated(enrichment_results_RNA)]


tss <- merge(tss, enrichment_results_RNA, by = "geneID", all.x = T)
tss_df <- tss[ !is.na(TSscore) , c("geneID", "tissue", "TSscore", "RNA_category" )] # subset only trusted values

unique(tss$RNA_category)
#tss_df[, category := "other"]
tss_df[RNA_category == "rna_hk", category := "house-keeping"] # this is how we can assign category
tss_df[TSscore >= 2.5, category := "tissue enriched"] # this is how we can assign category
tss_df[TSscore >= 4 & RNA_category == "rna_specific" , category := "tissue specific"] # this is how we can assign category


tss_df[ , dis_gene := ( geneID %in% unique(hpo_gene$geneID) ) ]
tss_df[is.na(category), category := "other"]


ggplot(tss_df, aes( x =  dis_gene, y=  TSscore, fill = dis_gene))+
  geom_boxplot() +
  facet_wrap( ~tissue )+
  theme_bw()

hist(tss_df$TSscore)


##############################
#### Analyse disease genes 
dg_tss <- tss_df[dis_gene == T ]
uniqueN(hpo_gene$geneID)
uniqueN(dg_tss$geneID)

dg_tss_gc <- dg_tss[, c("geneID", "category")]
dg_tss_gc[is.na(category), category := "other"]
dg_tss_gc <- dg_tss_gc[!duplicated(dg_tss_gc)]
dg_tss_gc <- dg_tss_gc[!is.na(category)] 

os <- dg_tss_gc[ , .N, by = category]

os[ , proportion := N/uniqueN(dg_tss_gc$geneID)]

##############################
#### annotate with hpo

# TODO use ontologyIndex function to get umrella_HPO_term
# head(pat_hpo)

head(hpo_gene)
dim(hpo_gene)
hpo_gene <- hpo_gene[geneID %in% unique(tss_df$geneID)]
hpo_gene <- hpo_gene[HPO_ID %in% hpo$id] # filter and remove ID that are not in ontology


hpo_ID_to_hpo_term <- hpo_gene[, c("HPO_ID")] 
hpo_ID_to_hpo_term <- hpo_ID_to_hpo_term[!duplicated(hpo_ID_to_hpo_term)]
hpo_ID_to_hpo_term$umbrella_hpo_term <- lapply(hpo_ID_to_hpo_term$HPO_ID, function(hpo_id) {
  get_term_property(ontology=hpo, property="ancestors", term=hpo_id, as_names=TRUE)[3]
})

hpo_ID_to_hpo_term <- hpo_ID_to_hpo_term %>% 
  mutate_at('umbrella_hpo_term', paste)
str(hpo_ID_to_hpo_term$umbrella_hpo_term)

hpo_gene_with_terms <- merge(hpo_gene, hpo_ID_to_hpo_term, by = "HPO_ID")

unique(hpo_gene_with_terms$umbrella_hpo_term)[!is.na(unique(hpo_gene_with_terms$umbrella_hpo_term))]

map_hpo_to_tissue <- read_excel("/Users/alexandraprimak/Desktop/hpo_to_tissue.xlsx" ,  
                                sheet = 1) %>% as.data.table

hpo_gene_with_terms <- merge(hpo_gene_with_terms, map_hpo_to_tissue, by = "umbrella_hpo_term", allow.cartesian=TRUE)
hpo_gene_with_terms <- hpo_gene_with_terms
tss_df_disease_gene <- tss_df[tss_df$dis_gene]


########### 

hpo_gene_with_terms$gene_tissue <- paste0(hpo_gene_with_terms$geneID, "_", hpo_gene_with_terms$tissue )

tss_df$gene_tissue <- paste0(tss_df$geneID, "_", tss_df$tissue )


tss_df[ , affected_tissue := gene_tissue %in% unique(hpo_gene_with_terms$gene_tissue )]
tss_df_disease_gene <- tss_df[tss_df$dis_gene] # matching

# Create new dataset
tis_dis <- tss_df_disease_gene[, 1:8]

# we need to remove tissue information 
unique(tis_dis$category)
tis_dis$category <- factor(tis_dis$category, 
                           levels = c("other", "house-keeping" , "tissue enriched", "tissue specific" ))

tis_dis <- tis_dis[order( category, decreasing = T)] # sort by category, assign importance by the level of factor in previous step 
head(tis_dis)


## One way to look on the data is like this
os <- tis_dis[, .N, by = .(category, affected_tissue)]
os[affected_tissue == T, total := nrow(tis_dis[affected_tissue == T])]
os[affected_tissue == F, total := nrow(tis_dis[affected_tissue == F])]
os[ , prop := N/ total] # 

p1 <- ggplot(os, aes(affected_tissue , prop, fill = category)) + 
  geom_bar(stat = 'identity') +
  coord_flip()+
  scale_y_continuous(labels=scales::percent) + 
  scale_fill_pander() + 
  labs(y = 'Proportion of (affected tissue) genes, showing expression pattern in any tissue', x = 'Affected tissue class') +
  theme_bw(base_size = 14)+
  theme(legend.position="top", 
        legend.title = element_blank())
p1


## another way to look at the same data 
vo <- table(tis_dis[,.(category, affected_tissue)]) %>% as.data.table()
vo[, prop := N/sum(N), by = category]
vo[, total := sum(N), by = category]
vo[, c("ci_l", "ci_h") := as.list(prop.test(N, total, correct = F)$conf.int), by = 1:nrow(vo)]
totals_vars_class <- vo[, sum(N), by = category]


p3 <- ggplot(vo[affected_tissue != F], aes(category, prop, fill = affected_tissue)) + geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_h), position = position_dodge(.9), width = .3) +
  scale_y_continuous(labels=scales::percent) + 
  scale_fill_pander() + 
  labs(y = 'Proportion of genes', x = 'Tissue specifisity class') +
  geom_text(aes(category, .18, label = V1, fill = NULL), data = totals_vars_class) +
  theme_bw(base_size = 14)
p3
# Indeed, tissue specific show enrichment, can we assign p-values here? 


## Two ways to calculate pvalues for categorial varaibles :
# First, Fisher exact test 
# Second - logistic regression ( glm( dis_gene ~ category , data = tis_dis , family = "binomial"))

## Do the enrichment via Fisher's test - short script 
tis_dis_df <- tis_dis
tis_dis_df$category <- as.character(tis_dis_df$category)
categories <- unique(tis_dis_df$category)

enr_fisher = data.table()
for(cat in categories){
  for(dc in unique(tis_dis_df$affected_tissue)){
    pt <- tis_dis_df
    pt <- pt[ , e_cat := category == cat]
    pt[, d_class := affected_tissue == dc]
    A = nrow(pt[d_class == T &  e_cat == T]) # disease gene in category  
    B = nrow(pt[d_class == T &  e_cat == F]) # disease gene NOT in category
    C = nrow(pt[d_class == F &  e_cat == T]) # NOT disease gene in category 
    D = nrow(pt[d_class == F &  e_cat == F]) # NOT disease gene NOT in category 
    
    risk_ratio = (A/(A+B))/(C/(C+D))
    odds_ratio = (A/B) / (C/D)
    
    mat <- table(pt[,.(d_class, e_cat)])
    ft <- fisher.test(matrix(c(A, B, C, D), ncol = 2))
    dt <- data.table(affected_tissue = dc, category = cat, fisher_pv = ft$p.value, 
                     RR = risk_ratio, OR = odds_ratio, log_odds = log(odds_ratio),  
                     ci_l = ft$conf.int[1], ci_h = ft$conf.int[2], 
                     log_ci_l = log(ft$conf.int[1]), log_ci_h = log(ft$conf.int[2]))
    enr_fisher <- rbind(enr_fisher, dt)
  }}

enr_fisher
setorderv(enr_fisher, c("affected_tissue", "fisher_pv"))
enr_fisher[ , signif := fisher_pv < 0.05]
enr_fisher[signif == T]


library(cowplot)

# Visualise enrichment results
p2 <- ggplot(enr_fisher[affected_tissue != F], aes(reorder(category, OR), log_odds)) + 
  geom_point() +
  ggtitle("Affected tissues")+ 
  geom_errorbar(aes(ymin = log_ci_l, ymax = log_ci_h), position = position_dodge(.9), width = 0) +
  geom_hline(yintercept = 0) + 
  labs(x = 'Gene category', y = 'Log Odds Ratio') + 
  theme_cowplot() + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank())
p2


library(patchwork)
p2 + p3 + plot_layout(widths = c(1, 2))

enrich <- function(data, outlier, features){
  estims = c()
  data <- as.data.frame(data)
  for(i in 1:length(features)){
    mod = summary(glm(data[, outlier]  ~ scale(data[, features[i]]), data = data, family = binomial))$coefficients
    estims.temp = data.frame(Cat = features[i], Estim = mod[2, 1], Std = mod[2, 2], Pval = mod[2, 4], stringsAsFactors = F)
    estims = rbind(estims, estims.temp)
  }
  estims = estims[order(estims$Estim, -estims$Pval), ]
  cat.names = as.character(estims$Cat)
  estims$Cat = factor(as.character(estims$Cat), levels = as.character(estims$Cat), labels = cat.names)
  return(estims)
}

enrich(tis_dis, "affected_tissue", c("TSscore"))


os <- tss_df[, .N, by = .(category, dis_gene)]
os[dis_gene == T, total := nrow(tss_df[dis_gene == T])]
os[dis_gene == F, total := nrow(tss_df[dis_gene == F])]
os[ , prop := N/ total]

p1 <- ggplot(os, aes( dis_gene , prop, fill = category)) + 
  geom_bar(stat = 'identity') +
  coord_flip()+
  scale_y_continuous(labels=scales::percent) + 
  scale_fill_pander() + 
  labs(y = 'Proportion of (disease) genes, showing expression pattern in any tissue', x = 'Tissue specifisity class') +
  theme_bw(base_size = 14)+
  theme(legend.position="top", 
        legend.title = element_blank())
p1
#We observe slight enrichment of tissue enriched expression in disease genes 
# now diference for other categories


## another way to look at the same data 
vo <- table(tss_df[,.(category, dis_gene)]) %>% as.data.table()
vo[, prop := N/sum(N), by = category]
vo[, total := sum(N), by = category]
vo[, c("ci_l", "ci_h") := as.list(prop.test(N, total, correct = F)$conf.int), by = 1:nrow(vo)]
totals_vars_class <- vo[, sum(N), by = category]


p3 <- ggplot(vo[dis_gene != F], aes(category, prop, fill = dis_gene)) + geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_h), position = position_dodge(.9), width = .3) +
  scale_y_continuous(labels=scales::percent) + 
  scale_fill_pander() + 
  labs(y = 'Proportion of genes', x = 'Tissue specifisity class') +
  geom_text(aes(category, .18, label = V1, fill = NULL), data = totals_vars_class) +
  theme_bw(base_size = 14)
# Indeed, tissue enriched show enrichment, can we assign p-values here? 



## Two ways to calculate pvalues for categorial varaibles :
# First, Fisher exact test 
# Second - logistic regression ( glm( dis_gene ~ category , data = tis_dis , family = "binomial"))

## Do the enrichment via Fisher's test - short script 
tis_dis_df <- tss_df
tis_dis_df$category <- as.character(tis_dis_df$category)
categories <- unique(tis_dis_df$category)

enr_fisher = data.table()
for(cat in categories){
  for(dc in unique(tis_dis_df$dis_gene)){
    pt <- tis_dis_df
    pt <- pt[ , e_cat := category == cat]
    pt[, d_class := dis_gene == dc]
    A = nrow(pt[d_class == T &  e_cat == T]) # disease gene in category  
    B = nrow(pt[d_class == T &  e_cat == F]) # disease gene NOT in category
    C = nrow(pt[d_class == F &  e_cat == T]) # NOT disease gene in category 
    D = nrow(pt[d_class == F &  e_cat == F]) # NOT disease gene NOT in category 
    
    risk_ratio = (A/(A+B))/(C/(C+D))
    odds_ratio = (A/B) / (C/D)
    
    mat <- table(pt[,.(d_class, e_cat)])
    ft <- fisher.test(matrix(c(A, B, C, D), ncol = 2))
    dt <- data.table(dis_gene = dc, category = cat, fisher_pv = ft$p.value, 
                     RR = risk_ratio, OR = odds_ratio, log_odds = log(odds_ratio),  
                     ci_l = ft$conf.int[1], ci_h = ft$conf.int[2], 
                     log_ci_l = log(ft$conf.int[1]), log_ci_h = log(ft$conf.int[2]))
    enr_fisher <- rbind(enr_fisher, dt)
  }}

enr_fisher
setorderv(enr_fisher, c("dis_gene", "fisher_pv"))
enr_fisher[ , signif := fisher_pv < 0.05]
enr_fisher[signif == T]

p2 <- ggplot(enr_fisher[dis_gene != F], aes(reorder(category, OR), log_odds)) + 
  geom_point() +
  ggtitle("Disease genes")+ 
  geom_errorbar(aes(ymin = log_ci_l, ymax = log_ci_h), position = position_dodge(.9), width = 0) +
  geom_hline(yintercept = 0) + 
  labs(x = 'Gene category', y = 'Log Odds Ratio') + 
  theme_cowplot() + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank())
p2

p1 + p2 + p3 + plot_layout(widths = c(c(1, 1, 1)))
p2 + p3 + plot_layout(widths = c(c(1, 1)))
# from this results, we see that disease genes likely will show some tissue enriched or specific expression!
# as expected, non- disease gene are enriched for category other, that is equal to: disease genes are depleted for category other!


# hmdata <- tis_dis[, c('geneID', 'TSscore')]
# 
# uhpoterms_gene <- hpo_gene_with_terms[, c("geneID", "tissue", "umbrella_hpo_term")]
# hmdata <- hmdata[geneID %in% unique(uhpoterms_gene$geneID)]
# uhpoterms_gene <- uhpoterms_gene[geneID %in% unique(hmdata$geneID)]
# 
# hmdata1 <- merge(hmdata, uhpoterms_gene, by="geneID", allow.cartesian=TRUE)
# hmdata1 <- hmdata1[!duplicated(hmdata1)]
# hmdata1[, meanTS := mean(TSscore, na.rm = T), by=c("umbrella_hpo_term", "tissue")]
# 
# hmdata1$geneID <- NULL
# hmdata1$TSscore <- NULL
# as.data.table(hmdata1)
# 
# hmdata1 <- hmdata1[!duplicated(hmdata1)]
# 
# hm <- dcast(hmdata1, umbrella_hpo_term ~ tissue)
# library(pheatmap)
# pheatmap(hm)
# 
# tss_dfX <- merge(tss_df, hpo_gene2, by = "geneID", allow.cartesian=TRUE) #all = T
# tss_dfX <- tss_dfX[!duplicated(tss_dfX)]

