library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(khroma)

model <- read_rds("./input/best_model.rds")
betas <- read_rds("./input/inspire/epigenomics.rds")
probes_450k <- read_rds("./input/450K_probes.rds")
betas <- betas[intersect(rownames(betas),probes_450k),]
betas <- betas[model$gene[2:nrow(model)],]
ic <- read_rds("./input/inspire/ic_scores.rds") %>% na.omit 
ic <- ic %>% filter(identAno%in%colnames(betas))
ic$identAno <- as.character(ic$identAno)
ic <- ic %>% arrange(IC)
betas <- betas[,ic$identAno]
cors <- apply(betas,1,function(x) cor(ic$IC,x, method = "s"))
cors <- cors %>% enframe %>% arrange(desc(value))
betas <- betas[cors$name,]

smooth_rainbow <- color("sunset")

col_fun <- colorRamp2(seq(0,1,0.01),smooth_rainbow(101, range = c(0, 1)))
col_ic <- colorRamp2(seq(0.5,1,0.0625),brewer.pal(9, "Greens"))
col_cor_ic <- colorRamp2(seq(-1,1,0.2),brewer.pal(11, "RdBu"))

ha = columnAnnotation(IC = ic$IC, col = list(IC = col_ic), show_legend = FALSE)
hb = rowAnnotation(`r IC` = cors$value, col = list(`r IC` = col_cor_ic), show_legend = FALSE)

ht <- Heatmap(betas,
              top_annotation = ha,
              left_annotation = hb,
              border = TRUE,
              row_title = "CpG",
              column_title = "Sample",
              cluster_columns = FALSE,
              row_dend_side = "right",
              row_names_side = "left",
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              cluster_rows = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              show_row_dend = FALSE,
              show_heatmap_legend = FALSE,
              show_column_dend = FALSE,
              heatmap_legend_param = list(direction = "vertical", title = expression("r"[s]), title_position = "topcenter", legend_height = unit(4, "cm")),
              col = col_fun) +
  rowAnnotation(label = anno_mark(padding = 3, at = c(1:5,87:91), labels = rownames(betas)[c(1:5,87:91)]))
a <- grid.grabExpr(draw(ht))

data <- read_rds("./input/best_model_predictions.rds")
data$sample <- as.numeric(data$sample)
ic <- read_rds("./input/inspire/ic_scores.rds") %>% na.omit
data <- data %>% left_join(ic[,c("identAno","age")], by = c("sample"="identAno"))
b <- ggplot(data, aes(ic, predicted_ic))+
  geom_point(fill = "#0077bb", shape = 21)+
  stat_correlation(label.x = 0.95, label.y = 0.05, method = "spearman", aes(label = paste0("r[s] == ",round(after_stat(rho),2))))+
  geom_smooth(method = "lm", color = "black")+
  theme_pubr(border = FALSE)+
  labs(x = "Intrinsic capacity", y = "DNAm intrinsic capacity")+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))

dnam <- read_rds("./input/inspire/epigenetic_clocks.rds")
combined <- dnam %>% column_to_rownames(var = "identAno")
pairwise <- reshape2::melt(combined[,c(1:3,5)] %>% as.matrix) %>% 
  purrr::set_names("identAno","variable","value") %>% 
  mutate(identAno = as.character(identAno)) %>% 
  left_join(dnam[,c("identAno","IC")])

pairwise %>% group_by(variable) %>% summarise(r = cor.test(IC,value, method = "s")$estimate, p = cor.test(IC,value, method = "s")$p.value)
d <- ggplot(pairwise, aes(IC, value))+
  geom_point(color = "grey",shape = 21)+
  geom_smooth(color = "black", method = "lm", se = FALSE)+
  stat_correlation(label.x = 0.5, label.y = 0.95, method = "spearman", aes(label = paste0("r[s] == ",round(after_stat(rho),2))))+
  theme_pubr(border = TRUE)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  facet_wrap(.~variable, nrow = 1, scales = "free_y")+
  labs(x = "DNAm intrinsic capacity", y = "DNAm age acceleration")+
  theme(legend.position = "none")

cor_mat <- cor(combined[,c(1:5)], use = "pairwise.complete.obs", method = "spearman")
col_fun <- colorRamp2(seq(-1,1,0.2),rev(brewer.pal(11, "RdBu")))
rownames(cor_mat) <- colnames(cor_mat) <- c("Horvath","Hannum","PhenoAge","IC","GrimAge")
ht <- Heatmap(cor_mat,
              rect_gp = gpar(col = "black", lwd = 1, lty = 2),
              border = TRUE,
              cell_fun = function(j, i, x, y, w, h, col) { 
                grid.text(round(cor_mat,2)[i, j], x, y)
              },
              cluster_columns = TRUE,
              row_dend_side = "right",
              row_names_side = "left",
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              heatmap_legend_param = list(direction = "vertical", title = expression("r"[s]), title_position = "topcenter", legend_height = unit(4, "cm")),
              col = col_fun)
c <- grid.grabExpr(draw(ht))

library(ggvenn)
library(gridExtra)

model <- read_rds("./input/best_model.rds")
coefIC <- model$gene[2:nrow(model)]
coefHannum <- read_rds("./input/clocks/coefHannum.rds")
coefHorvath <- read_rds("./input/clocks/coefHorvath.rds")
coefPhenoAge <- read_rds("./input/clocks/coefLevine.rds")
coefGrimAge2 <- read_rds("./input/clocks/coefGrimAge2.rds")
coefPhenoAge <- coefLevine$CpGmarker
coefPhenoAge <- unique(coefPhenoAge[grepl("cg",coefPhenoAge)])
coefHannum <- coefHannum$CpGmarker
coefHannum <- unique(coefHannum[grepl("cg",coefHannum)])
coefHorvath <- coefHorvath$CpGmarker
coefHorvath <- unique(coefHorvath[grepl("cg",coefHorvath)])

cpgs <- list(Hannum = coefHannum,
             Horvath = coefHorvath,
             PhenoAge = coefPhenoAge,
             GrimAge = coefGrimAge2,
             IC = coefIC)

p1 <- ggvenn(cpgs[c(5,1)],
             show_percentage = FALSE,
             stroke_color = c(rep("#004488",100),rep("#DDAA33",100)),
             fill_color = rep("white",2),
             fill_alpha = 0.5,
             stroke_size = 2, set_name_size = 4) +
  coord_cartesian(clip="off")

p2 <- ggvenn(cpgs[c(5,2)],
             show_percentage = FALSE,
             stroke_color = c(rep("#004488",100),rep("#DDAA33",100)),
             fill_color = rep("white",2),
             fill_alpha = 0.5,
             stroke_size = 2, set_name_size = 4) +
  coord_cartesian(clip="off")

p3 <- ggvenn(cpgs[c(5,3)],
             show_percentage = FALSE,
             stroke_color = c(rep("#004488",100),rep("#DDAA33",100)),
             fill_color = rep("white",2),
             fill_alpha = 0.5,
             stroke_size = 2, set_name_size = 4) +
  coord_cartesian(clip="off")

p4 <- ggvenn(cpgs[c(5,4)],
             show_percentage = FALSE,
             stroke_color = c(rep("#004488",100),rep("#DDAA33",100)),
             fill_color = rep("white",2),
             fill_alpha = 0.5,
             stroke_size = 2, set_name_size = 4) +
  coord_cartesian(clip="off")

blank <- grid.rect(gp=gpar(col="white"))
e <- grid.arrange(
  blank,blank,blank,blank,blank,blank,blank,blank,blank,
  blank,p1,blank,p2,blank,p3,blank,p4,blank, nrow = 2, widths = c(0.1,1,0.2,1,0.2,1,0.2,1,0.1), heights = c(0.1,1))

pdf(file = "./output/figures/figure2.pdf",width=12,height=8.70)
plot_grid(
  plot_grid(a, b, c, nrow = 1, labels = c("a","b","c")),
  plot_grid(d, NULL, nrow = 1, labels = c("d"), rel_widths = c(1,0.01)),
  plot_grid(NULL,e, NULL, nrow = 1, labels = c("e"), rel_widths = c(0.04,1,0.02)),
  nrow = 3,
  rel_heights = c(1,0.8,0.45),
  label_size = 12,
  align = "v"
)
dev.off()
