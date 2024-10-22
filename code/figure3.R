library(ggrepel)
library(CellPlot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(fgsea)
library(dplyr)
library(impute)

raw <- read_rds("./input/fhs/epigenomics.rds") 
raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
data_adjusted <- data %>% t()
predictions <- apply(data_adjusted, 2, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8, sex, date8) %>% na.omit %>% set_names("sample","age","sex","date")
ages$sample <- as.character(ages$sample)
ages$sex <- ifelse(ages$sex==1,"male","female")
predictions <- predictions %>% left_join(ages) %>% na.omit
predictions$ic_acc <- resid(lm(ic~age, data = predictions))

meta <- read_tsv("./input/fhs/phe000002.v10_release_manifest.txt", comment = "#")
meta <- meta[,1:2] %>% set_names("file","sample")
meta$sample <- as.character(meta$sample)
predictions <- predictions %>% left_join(meta) %>% na.omit

exp <- read_tsv("./input/fhs/FinalFile_Gene_OFF_2446_Adjusted_c1.txt.gz")
exp_long <- reshape2::melt(exp %>% column_to_rownames(var = "transcript_cluster_id") %>% as.matrix)
exp_long$Var1 <- as.character(exp_long$Var1)
exp_long$Var2 <- as.character(exp_long$Var2)
exp_long <- exp_long %>% filter(Var2%in%unique(predictions$file))

annotation <- read_tsv("input/fhs/GPL5175.txt.gz", skip = 14)
annotation$gene <- sapply(annotation$gene_assignment, function(x) strsplit(x, split = "\\ // ")[[1]][2])
annotation <- annotation %>% filter(ID%in%unique(exp_long$Var1))

exp_mat <- exp_long %>% left_join(annotation[,c("ID","gene")], by = c("Var1"="ID")) %>% na.omit
exp_mat <- exp_mat %>% group_by(gene,Var2) %>% summarise(value = mean(value))
exp_mat <- reshape2::acast(exp_mat, gene~Var2, value.var = "value")

predictions <- predictions %>% filter(file%in%colnames(exp_mat))
exp_mat <- exp_mat[,predictions$file]

deg_ic <- apply(exp_mat, 1, function(x) summary(lm(predictions$ic_acc~x))$coefficients[2,]) %>% t() %>% data.frame %>% set_names("estimate","se","t","p") %>% rownames_to_column(var = "gene")
deg_ic$fdr <- p.adjust(deg_ic$p)
deg_ic$gene[deg_ic$gene=="TRA@"] <- "TCRA"
deg_ic %>% filter(fdr<0.05) %>% nrow

a <- ggplot(deg_ic, aes(estimate, -log10(p)))+
  geom_point(fill = "grey", shape = 21, color = "grey")+
  geom_point(data = deg_ic %>% filter(fdr<0.05&estimate>0), fill = "#d6604d", shape = 21)+
  geom_point(data = deg_ic %>% filter(fdr<0.05&estimate<0), fill = "#4393c3", shape = 21)+
  theme_pubr()+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_text_repel(data = deg_ic %>% filter(fdr<0.05&estimate<0) %>% arrange(p) %>% slice_head(n = 20), aes(label = gene), size = 3, max.overlaps = 50)+
  geom_text_repel(data = deg_ic %>% filter(fdr<0.05&estimate>0) %>% arrange(p) %>% slice_head(n = 20), aes(label = gene), size = 3, max.overlaps = 50)+
  labs(x = "Log2 Fold Change", y = "-Log10 (P-value)")+
  theme(plot.background=element_blank(),strip.background = element_blank(), strip.text = element_text(size = 14))

b <- ggdraw() + draw_image(magick::image_read_pdf("./input/3b.pdf", pages = 1)) 

geneset <- read_rds("./input/hoa_genesets.rds")
signature <- deg_ic
signature$score <- -log10(signature$p)
signature <- setNames(signature$score, signature$gene)
results <- fgseaSimple(pathways = geneset, stats = signature, scoreType = "pos", nperm = 1000000)
results$score <- -log10(results$pval)
results$pathway <- str_to_sentence(results$pathway)
results$pathway <- gsub(" ","\n",results$pathway)
results <- results %>% mutate(pathway = factor(pathway, levels = unique(pathway)))

c <- ggplot(results, aes(x = pathway, y = score, fill = score)) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(clip = "off") +  
  geom_hline(yintercept = max(results$score) * 1.05, color = "black", linetype = "solid") + 
  geom_segment(data = results, aes(x = as.numeric(pathway), 
                                   xend = as.numeric(pathway), 
                                   y = 0, 
                                   yend = max(results$score) * 1.1),
               linetype = "dotted", color = "black") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_blank(),  
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",  
    legend.justification = "center",  
    legend.title = element_text(hjust = 0.5, size = 12),  
    plot.margin = unit(c(0, 0, 0, 0), "cm") 
  ) +
  scale_y_continuous(limits = c(0, max(results$score) * 1.2)) +  
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd")) + 
  labs(x = "", y = "", fill = "-log10(p-value)")

#figure d
cellCountPredictorHorvath=function(dat0, datCOEF,imputeValues=TRUE) {
  datout=data.frame(matrix(NA,nrow=dim(dat0)[[2]]-1,ncol=dim(datCOEF)[[2]]-1 ))
  match1=match(datCOEF[-1,1, drop = TRUE],dat0[,1] )
  if (   sum(!is.na(match1))==0 ) stop("Input error. The first column of dat0 does not contain CpG identifiers (cg numbers).")   
  dat1=dat0[match1,]
  row.names1=as.character(dat1[,1])
  dat1=dat1[,-1]
  if (imputeValues ){dat1=impute.knn(data=as.matrix(dat1) ,k = 10)[[1]]}
  for (i in 1:dim(dat1)[[2]] ){ 
    for (j in 2:dim(as.matrix(datCOEF))[[2]] ){
      datout[i,j-1]=sum(dat1[,i]* datCOEF[-1,j],na.rm=TRUE)+ datCOEF[1,j]}
  }
  colnames(datout)=colnames(datCOEF)[-1]
  rownames(datout)=colnames(dat0)[-1]
  datout=data.frame(SampleID= colnames(dat0)[-1],datout)
  datout
} 

datCellCountPredictorsHorvathSubset <- read_rds("./input/cell_count_predictor_horvath.rds")
raw <- read_rds("./input/fhs/epigenomics.rds") 
dat1 <- raw %>% data.frame(check.names = FALSE) %>% rownames_to_column(var = "Probe")
datCellCountsHorvath=cellCountPredictorHorvath(dat1,datCOEF=datCellCountPredictorsHorvathSubset,imputeValues=TRUE) 
colnames(datCellCountsHorvath)[1] <- "sample"
cell_counts <- methylclock::meffilEstimateCellCountsFromBetas(beta = raw, cellTypeReference = "blood gse35069", verbose = TRUE)
cell_counts <- cell_counts %>% data.frame %>% rownames_to_column(var = "sample")
cell_counts <- cell_counts %>% left_join(datCellCountsHorvath)
colnames(cell_counts) <- c("sample","B","CD4+", "CD8+","Granulocyte","Monocyte","NK","Plasmablast", "CD8+ Exhausted","CD8+ Naive")

raw <- raw %>% t()
model <- read_rds("./input/best_model.rds")
data <- raw[,model$gene[2:nrow(model)]] 
data_adjusted <- data %>% t()
predictions <- apply(data_adjusted, 2, function(x) sum(x*model$coefficient[2:nrow(model)])+model$coefficient[1]) %>% 
  enframe %>%
  set_names("sample","ic") %>% 
  data.frame
predictions <- predictions %>% left_join(cell_counts)
ages <- read_tsv("./input/fhs/phs000007.v33.pht003099.v8.p14.c1.vr_dates_2019_a_1175s.HMB-IRB-MDS.txt", comment = "#")
ages <- ages %>% dplyr::select(shareid, age8) %>% na.omit %>% set_names("sample","age")
ages$sample <- as.character(ages$sample)
predictions <- predictions %>% left_join(ages) %>% na.omit 

age_changes <- reshape2::melt(predictions[,c(1,3:11)]) %>% left_join(predictions[,c(1,12)])
age_changes %>% group_by(variable) %>%  summarise(r = cor.test(age,value)$estimate, p = cor.test(age,value)$p.value)

p1 <- ggplot(age_changes, aes(age,value))+
  geom_point(shape = 21, fill = "grey", alpha = 0.1)+
  geom_smooth(method = "lm", color = "red")+
  stat_correlation(label.x = 0.05, label.y = 0.95, boot.R = 0, method = "spearman", exact = TRUE, color = "black",aes(label = paste0("r[s] == ",round(after_stat(rho),2))))+
  stat_correlation(label.x = 0.05, label.y = 0.85, boot.R	= 0, method = "spearman", exact = TRUE, color = "black", aes(label = paste0("p == ",after_stat(sprintf("%.2e", p.value)))))+
  labs(x = "Age", y = "Cell count")+
  facet_wrap(.~variable, nrow = 1, scales = "free_y")+
  theme_pubr(border = TRUE)

predictions$ic <- resid(lm(ic~age, data = predictions))
predictions$B <- resid(lm(B~age, data = predictions))
predictions$`CD4+` <- resid(lm(`CD4+`~age, data = predictions))
predictions$`CD8+` <- resid(lm(`CD8+`~age, data = predictions))
predictions$Granulocyte <- resid(lm(Granulocyte~age, data = predictions))
predictions$Monocyte <- resid(lm(Monocyte~age, data = predictions))
predictions$NK <- resid(lm(NK~age, data = predictions))
predictions$Plasmablast <- resid(lm(Plasmablast~age, data = predictions))
predictions$`CD8+ Exhausted` <- resid(lm(`CD8+ Exhausted`~age, data = predictions))
predictions$`CD8+ Naive` <- resid(lm(`CD8+ Naive`~age, data = predictions))

cell_changes <- reshape2::melt(predictions[,c(1,3:11)]) %>% left_join(predictions[,c(1,2)])
cell_changes %>% group_by(variable) %>%  summarise(r = cor.test(ic,value)$estimate, p = cor.test(ic,value)$p.value)
p2 <- ggplot(cell_changes, aes(ic,value))+
  geom_point(shape = 21, fill = "grey", alpha = 0.1)+
  geom_smooth(method = "lm", color = "red")+
  stat_correlation(label.x = 0.05, label.y = 0.95, boot.R = 0, method = "spearman", exact = TRUE, color = "black",aes(label = paste0("r[s] == ",round(after_stat(rho),2))))+
  stat_correlation(label.x = 0.05, label.y = 0.85, boot.R	= 0, method = "spearman", exact = TRUE, color = "black", aes(label = paste0("p == ",after_stat(sprintf("%.2e", p.value)))))+
  labs(x = "DNAm IC (age-adjusted)", y = "Cell count (age-adjusted)")+
  facet_wrap(.~variable, nrow = 1, scales = "free_y")+
  theme_pubr(border = TRUE)

d <- grid.arrange(p1,p2, nrow = 2)

#combine figures
plot_grid(
  plot_grid(a, b, c, nrow = 1, labels = c("a","b","c"), rel_widths = c(0.8,1,1)),
  plot_grid(d, nrow = 1, labels = c("d")),
  nrow = 2,
  label_size = 12,
  rel_heights = c(1,1),
  align = "v"
)
