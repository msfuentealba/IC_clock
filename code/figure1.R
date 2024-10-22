library(tidyverse)
library(ggpubr)
library(chngpt)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

age <- read_delim("./input/inspire/var_calc.csv")[,c("identAno","ageP1M0")] %>% set_names("identAno","age")
sex <- read_delim("./input/inspire/sex_dateNaissP.csv")[,c("identAno","sexeP")] %>% set_names("identAno","sex")
sex$sex <- ifelse(sex$sex=="2","Female","Male")
females <- sex$identAno[sex$sex=="Female"]

#mmse
mmse <- read_delim("./input/inspire/mmse.csv")
mmse <- mmse[,c(1,9:ncol(mmse))]
mmse[mmse=="D"] <- NA
mmse <- mmse %>% mutate_all(as.numeric)

#sbbp
sppb <- read_delim("./input/inspire/var_calc.csv")
sppb <- sppb[,c("identAno","SPPB1scM0","SPPB2scM0","SPPB3scM0")]

#psychological
qsp <- read_delim("./input/inspire/trbSens_QSP9.csv")
qsp <- qsp[,c(1,21:29)] 

#hand grip strength
hand <- read_delim("./input/inspire/Fried.csv")
hand <- hand[,c("identAno","FriedCpr1M0","FriedCpr2M0","FriedCpr3M0")]
hand <- hand[rowSums(is.na(hand[,2:ncol(hand)]))!=3,]

#sensory
eye_chart <- read_delim("./input/inspire/trbSens_QSP9.csv")
eye_chart <- eye_chart[,c("identAno","TBVIS1loinODM0","TBVIS1loinOGM0","TBVIS2loinODM0","TBVIS2loinOGM0","TBVIS3loinODM0","TBVIS3loinOGM0","TBVIS4presODM0","TBVIS4presOGM0")]
whisper <- read_delim("./input/inspire/st.csv") %>% filter(VISITNUM==0)
whisper <- whisper[,c("identAno","STAUDD","STAUDG")]

#calculate sum
cognitive <- tibble(identAno = mmse$identAno, MMSE = rowSums(mmse[,2:ncol(mmse)]))
locomotion <- tibble(identAno = sppb$identAno, SPPB = rowSums(sppb[,2:ncol(sppb)]))
psychology <- tibble(identAno = qsp$identAno, QSP = rowSums(qsp[,2:ncol(qsp)])) 
vitality <- tibble(identAno = hand$identAno, GRIPSTRENGTH = apply(hand[,2:ncol(hand)], 1, function(x) max(x, na.rm = TRUE))) 
sensory <- tibble(identAno = eye_chart$identAno, EYECHART = rowSums(eye_chart[,2:ncol(eye_chart)])) %>%
  left_join(tibble(identAno = whisper$identAno, WHISPER =  rowSums(whisper[,2:ncol(whisper)])))

#scale 0-1
cognitive_scaled <- cognitive %>% mutate(MMSE = (MMSE-0)/(30-0))
locomotion_scaled <- locomotion %>% mutate(SPPB = (SPPB-0)/(12-0))
psychology_scaled <- psychology %>% mutate(QSP = (QSP-27)/(0-27))
max_grip_female <- max(vitality$GRIPSTRENGTH[vitality$identAno%in%females], na.rm = TRUE)
max_grid_male <- max(vitality$GRIPSTRENGTH[!vitality$identAno%in%females], na.rm = TRUE)
vitality_scaled <- vitality %>% mutate(GRIPSTRENGTH = ifelse(identAno%in%females,(GRIPSTRENGTH-0)/(max_grip_female-0),(GRIPSTRENGTH-0)/(max_grid_male-0)))
sensory_scaled <- sensory %>% mutate(EYECHART = (EYECHART-0)/(8-0), WHISPER = (WHISPER-0)/(2-0))

#z-score 
ic <- age %>% left_join(sex) %>% left_join(cognitive_scaled) %>% left_join(locomotion_scaled) %>% left_join(psychology_scaled) %>% left_join(vitality_scaled) %>% left_join(sensory_scaled)
ic <- cbind(ic[,1:3],apply(ic[,4:9], 2, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)))
ic <- ic %>% mutate(identAno = identAno, age = age, sex = sex, Cognition = MMSE, Locomotion = SPPB, Psychological = QSP, Vitality = GRIPSTRENGTH, Sensory = (EYECHART+WHISPER)/2)
ic <- ic[,c("identAno","age","sex","Cognition", "Locomotion", "Psychological", "Vitality", "Sensory")]
ic$IC <- rowMeans(ic[,c("Cognition", "Locomotion", "Psychological", "Vitality", "Sensory")])
ic <- cbind(ic[,1:3],apply(ic[,4:ncol(ic)],2, function(x) (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))))
write_rds(ic, file = "./input/inspire/ic_scores.rds")

#reorder data
all_domains <- ic %>% na.omit
table(all_domains$sex)
ic_age <- reshape2::melt(ic[,c(1,4:9)] %>% column_to_rownames(var = "identAno") %>% as.matrix) %>% set_names("identAno","domain","score") %>% left_join(ic[,c("identAno","age","sex")]) %>% na.omit
ic_age$domain <- factor(ic_age$domain, levels = c("Cognition", 
                                                  "Locomotion", 
                                                  "Psychological", 
                                                  "Sensory",
                                                  "Vitality",
                                                  "IC"))
ptable <- ic_age %>% group_by(domain) %>% summarise(r = cor.test(age,score, method = "s")$estimate, p = cor.test(age,score, method = "s")$p.value) %>% mutate(fdr = p.adjust(p))

#figure 1b
b <- ggplot(ic_age, aes(age, score, fill = domain, color = domain))+
  geom_point(shape = 21, alpha = 0.2, size = 1.5)+
  geom_smooth(method = "loess")+
  stat_correlation(label.x = 0.05, label.y = 0.15, boot.R = 0, method = "spearman", exact = TRUE, color = "black",aes(label = paste0("r[s] == ",round(after_stat(rho),2))))+
  stat_correlation(label.x = 0.05, label.y = 0.05, boot.R	= 0, method = "spearman", exact = TRUE, color = "black", aes(label = paste0("p == ",after_stat(sprintf("%.2e", p.value)))))+
  theme_pubr(border = TRUE)+
  facet_wrap(.~domain, nrow = 2, scales = "free_x")+
  scale_fill_manual(values = c("#692e35","#182336","#5a8898","#c26236","#e9a04f","grey40"))+
  scale_color_manual(values = c("#692e35","#182336","#5a8898","#c26236","#e9a04f","grey40"))+
  labs(x = "Age", y = "Normalized scores", color = "Sex")+
  theme(legend.position = "none")

#figure 1c
c <- ggplot(ic_age %>% left_join(ic_age %>% group_by(sex,domain) %>% summarise(mean = mean(score))), aes(sex, score, color = sex))+
  geom_jitter(width = 0.2)+
  facet_wrap(.~domain)+
  geom_hline(aes(yintercept = mean), size = 1, color = "black")+
  geom_hline(aes(yintercept = mean, color = sex), size = 0.5)+
  stat_compare_means(comparisons = list(c("Female","Male")), vjust = -0.25)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0, 0.2)))+
  theme_pubr(border = TRUE)+
  scale_color_manual(values = c("#c2a5cf","#acd39e"))+
  theme(legend.position = "top", axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  labs(color = "Sex", x = "Sex", y = "Normalized score")

#figure 1d
breakpoint <- function(data){
  model <- chngptm(formula.1=score~1, formula.2=~age, data, type="M12c", family="gaussian", ci.bootstrap.size=100)
  line <- tibble(line = sapply(ic_age$age, function(x) (model$coefficients[1]) + (model$coefficients[2]*x) + (model$coefficients[3]*(ifelse(x>model$coefficients[4],x-model$coefficients[4],0))^2)) %>% as.numeric(),
                 age = ic_age$age,
                 coef = model$coefficients[4],
                 identAno = ic_age$identAno)
}
breakpoints <- ic_age %>% group_by(domain,sex) %>% nest() %>% mutate(breakpoint = map(data, breakpoint)) %>% dplyr::select(domain,sex,breakpoint) %>% unnest(breakpoint)
breakpoints <- breakpoints %>% dplyr::select(domain,sex,coef) %>% unique
breakpoints$domain <- factor(breakpoints$domain, levels = breakpoints %>% group_by(domain) %>% summarise(mean = mean(coef)) %>% arrange(mean) %>% pull(domain) %>% as.character() %>% rev)
breakpoints$nudge <- ifelse(breakpoints$sex=="Male",0.17,-0.17)

d <- ggplot(breakpoints, aes(coef, domain)) +
  geom_errorbarh(aes(xmin = 0, xmax = coef, color = sex),height = 0,position = position_dodge(width = .7)) +
  geom_point(aes(color = sex), size = 6, position = position_dodge(width = .7))+
  geom_text(aes(label = coef), parse = TRUE, position = position_nudge(y = breakpoints$nudge))+
  scale_color_manual(values = c("#c2a5cf","#acd39e"))+
  labs(x = "Age at decline", y = "Domain", color = "Sex")+
  theme_pubr()

#figure 1e
correlations <- cor(as.matrix(ic[,c("IC","Cognition","Locomotion","Vitality","Psychological","Sensory")]%>%na.omit), method = "s")
cor.test(ic$Locomotion,ic$Cognition, method = "s")$p.value
cor.test(ic$Psychological,ic$Sensory, method = "s")$p.value
col_fun <- colorRamp2(seq(-1,1,0.2),rev(brewer.pal(11, "RdBu")))
rownames(correlations) <- colnames(correlations) <- c("IC","Cognition","Locomotion","Vitality","Psychological","Sensory")
ht <- Heatmap(correlations,
        rect_gp = gpar(col = "black", lwd = 1, lty = 2),
        border = TRUE,
        cell_fun = function(j, i, x, y, w, h, col) { 
          grid.text(round(correlations,2)[i, j], x, y)
        },
        cluster_columns = TRUE,
        row_dend_side = "right",
        row_names_side = "left",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        cluster_rows = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        show_column_dend = TRUE,
        show_row_dend = FALSE,
        heatmap_legend_param = list(direction = "vertical", title = expression("r"[s]), title_position = "topcenter", legend_height = unit(4, "cm")),
        col = col_fun)
e <- grid.grabExpr(draw(ht))
  
#load figure 1a
a <- ggdraw() + draw_image(magick::image_read_pdf("./input/1a.pdf", pages = 1)) 

#combine figures
pdf(file = "./output/figures/figure1.pdf",width=13.40,height=10.04)
plot_grid(
  plot_grid(a, b, nrow = 1, labels = c("a","b"), rel_widths = c(1,2)),
  plot_grid(c,d,e,nrow = 1, labels = c("c","d","e")),
  nrow = 2,
  rel_heights = c(1.3,1),
  label_size = 12,
  align = "v"
)
dev.off()
