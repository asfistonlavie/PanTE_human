#### Characterization of Transposable Elements Variations in the Human Pangenome
#### Shadi Shahatit, Master's Thesis, ISEM, UM 2023-2024
# Libraries ---------------------------------------------------------------

library("qqman")
library("dplyr")
library("ggplot2")
library("ggpattern")
library("diveRsity")
library("rehh")
library("vcfR")
library("rehh.data")
library("writexl")
library("data.table")
library("RColorBrewer")
library("tidyverse")
library("forcats")
library("wesanderson")
library("stringr")
library("tidyr")
library("readr")
library("splitstackshape")
library("ggridges")
library("karyoploteR")
library("BiocManager")
library("GenomicRanges")
library("cowplot")
library("scales")
library("palmerpenguins")
library("ggbeeswarm")
library("systemfonts")
library("ggforce")
library("ggpubr")
library("rtracklayer")
library("ggcorrplot")
library("scMethrix")
library("BRGenomics")
library("minpack.lm")
library("sjPlot")
library("DHARMa")
library("RRPP")
library("devtools")
library("flexplot")
library("lme4")
library("betareg")
library("nnet")
library("MASS")
library("ordinal")
library("brant")
library("ggeffects")
library("lvplot")
library("ggthemes")
library("caret")
library("patchwork")
library("ggrepel")
library("MuMIn")
library("lmtest")
library("fitdistrplus")
library("VennDiagram")
library("viridis")
library("grid")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenometriCorr")
library("rGREAT")
library("regioneR")
library("plyranges")

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S3_Project"

# Figure 1 ----------------------------------------------------------------



# Define AFR, outofAFR, global allele freq, and TE annotation ----------------------------

## outofAFR freq

outofAFR_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/outofAFR/freqoutofAFR.frq"), header = T,
                            sep = " ",
                            dec = ".",
                            fill = T) 
colnames(outofAFR_freq)= c("V1")
outofAFR_freq$V1 <- gsub(":", "\t", outofAFR_freq$V1)
outofAFR_freq_2 <- str_split(outofAFR_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)
outofAFR_freq_3 <- outofAFR_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(outofAFR_freq_3) = c("CHR","POS","REF.o","FREQ_REF.o","ALT.o","FREQ_ALT.o")
outofAFR_freq_3[,6] <- as.numeric(unlist(outofAFR_freq_3[,6]))
outofAFR_freq_3[,4] <- as.numeric(unlist(outofAFR_freq_3[,4]))
outofAFR_freq_3 <- data.frame(append(outofAFR_freq_3, c(superpop.o="outofAFR"), after=6))

## AFR freq

AFR_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/AFR/freqAFR.frq"), header = T,
                       sep = " ",
                       dec = ".",
                       fill = T) 
colnames(AFR_freq)= c("V1")
AFR_freq$V1 <- gsub(":", "\t", AFR_freq$V1)
AFR_freq_2 <- str_split(AFR_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)
AFR_freq_3 <- AFR_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(AFR_freq_3) = c("CHR","POS","REF","FREQ_REF","ALT","FREQ_ALT")
AFR_freq_3[,6] <- as.numeric(unlist(AFR_freq_3[,6]))
AFR_freq_3[,4] <- as.numeric(unlist(AFR_freq_3[,4]))
AFR_freq_3 <- data.frame(append(AFR_freq_3, c(superpop="AFR"), after=6))

## global freq

super_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/superpop/freqSuper.frq"), header = T,
                         sep = " ",
                         dec = ".",
                         fill = T) 

colnames(super_freq)= c("V1")

super_freq$V1 <- gsub(":", "\t", super_freq$V1)

super_freq_2 <- str_split(super_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)

super_freq_3 <- super_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(super_freq_3) = c("CHR","POS","REF","FREQ_REF","ALT","FREQ_ALT")

super_freq_3[,6] <- as.numeric(unlist(super_freq_3[,6]))
super_freq_3[,4] <- as.numeric(unlist(super_freq_3[,4]))

super_freq_3[,2] <- as.numeric(unlist(super_freq_3[,2]))
super_freq_7 <- super_freq_3 %>% mutate(end=POS+nchar(ALT)-1)
super_freq_7 <- super_freq_7 %>% dplyr::select(c(1,2,6,7))

## TE_anno

TE_anno <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/RTE_norm_mm90_1alt.vcf"), header = F,
                      sep = " ",
                      dec = ".")
colnames(TE_anno) = c("CHR","POS","svid","ALT","ANNO")
TE_anno$ANNO <- gsub(";", "", TE_anno$ANNO)
TE_anno$ANNO <- gsub("=", "", TE_anno$ANNO)
TE_anno$ANNO <- gsub("rC", "\t", TE_anno$ANNO)
TE_anno[,6:7] <- stringr::str_split_fixed(TE_anno$ANNO, "\t", 2)
TE_anno <- TE_anno %>% dplyr::select(c(1,2,4,7))
TE_anno[,2] <- as.numeric(unlist(TE_anno[,2]))
TE_anno <- TE_anno[order(TE_anno[,1],TE_anno[,2]),]
colnames(TE_anno) = c("CHR","POS","ALT","ANNO")
TE_anno <- TE_anno %>% mutate(END=POS+nchar(ALT)-1)
super_freq_anno <- merge(super_freq_3,TE_anno)

TE_anno_Gr <- makeGRangesFromDataFrame(df=TE_anno,
                                       keep.extra.columns=T,
                                       ignore.strand=T,
                                       seqnames.field="CHR",
                                       start.field="POS",
                                       end.field="end",
                                       starts.in.df.are.0based=FALSE)

# TE size and Karyoplote ---------------------------------------------------

## TE size

chr_length <- read.table(file = file.path(sys_dir,"PanTE_human/karyoplote_and_TE_size/chr_lenght_hg28"), header = FALSE,
                         sep = "",
                         dec = ".")
super_freq_3$TE_size <- nchar(super_freq_3$ALT)
colnames(chr_length) <- c("CHR","chr_length")
chr_length$CHR <- gsub("chr01", "chr1", chr_length$CHR)
chr_length$CHR <- gsub("chr02", "chr2", chr_length$CHR)
chr_length$CHR <- gsub("chr03", "chr3", chr_length$CHR)
chr_length$CHR <- gsub("chr04", "chr4", chr_length$CHR)
chr_length$CHR <- gsub("chr05", "chr5", chr_length$CHR)
chr_length$CHR <- gsub("chr06", "chr6", chr_length$CHR)
chr_length$CHR <- gsub("chr07", "chr7", chr_length$CHR)
chr_length$CHR <- gsub("chr08", "chr8", chr_length$CHR)
chr_length$CHR <- gsub("chr09", "chr9", chr_length$CHR)
super_freq_3_chr_length <- merge(super_freq_3,chr_length) 

size_length <- super_freq_3_chr_length %>% dplyr::select(c("CHR","POS","chr_length","TE_size"))

size_length$TE_percen <- (size_length$TE_size/size_length$chr_length)*100
(sum(size_length$TE_percen))

size_length <- size_length %>% filter(TE_size >= 100)

TE_size_sum <- sum(size_length$TE_size)
chr_length_sum <- sum(unique(size_length$chr_length))
TE_genome_percentage <- ((TE_size_sum)/(chr_length_sum))*100
median(size_length$TE_size)
mean(size_length$TE_size)

size_length_counted <- size_length %>% dplyr::count(CHR,chr_length)
cor <- cor(size_length_counted$chr_length,size_length_counted$n)
ggplot(data=size_length_counted, mapping=aes(x=chr_length, y=n)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE count")+
  theme_classic() +
  annotate("text", x=min(size_length_counted$chr_length), y=max(size_length_counted$n), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEcountchrlength.png"), plot=ggplot2::last_plot())

TE_size_chr <- size_length %>%
  group_by(CHR) %>%
  summarise(total_TE_size = sum(TE_size),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()
cor <- cor(TE_size_chr$mean_chr_length,TE_size_chr$total_TE_size)
ggplot(data=TE_size_chr, mapping=aes(x=mean_chr_length, y=total_TE_size)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE size")+
  theme_classic() +
  annotate("text", x=min(TE_size_chr$mean_chr_length), y=max(TE_size_chr$total_TE_size), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsizechrlength.png"), plot=ggplot2::last_plot())

TE_percen_chr <- size_length %>%
  group_by(CHR) %>%
  summarise(total_TE_percen = sum(TE_percen),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()
cor <- cor(TE_percen_chr$mean_chr_length,TE_percen_chr$total_TE_percen)
ggplot(data=TE_percen_chr, mapping=aes(x=mean_chr_length, y=total_TE_percen)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE percentage")+
  theme_classic() +
  annotate("text", x=min(TE_percen_chr$mean_chr_length), y=max(TE_percen_chr$total_TE_percen), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsizepercenchrlength.png"), plot=ggplot2::last_plot())

kruskal_test <- kruskal.test(total_TE_percen ~ CHR, data = TE_percen_chr)
print(kruskal_test)

size_length$CHR <- gsub("chr1", "chr01", size_length$CHR)
size_length$CHR <- gsub("chr2", "chr02", size_length$CHR)
size_length$CHR <- gsub("chr3", "chr03", size_length$CHR)
size_length$CHR <- gsub("chr4", "chr04", size_length$CHR)
size_length$CHR <- gsub("chr5", "chr05", size_length$CHR)
size_length$CHR <- gsub("chr6", "chr06", size_length$CHR)
size_length$CHR <- gsub("chr7", "chr07", size_length$CHR)
size_length$CHR <- gsub("chr8", "chr08", size_length$CHR)
size_length$CHR <- gsub("chr9", "chr09", size_length$CHR)

# hist(size_length$TE_size, 
#      col = "skyblue",             
#      border = "black",            
#      main = "Histogram of Size",  
#      xlab = "Size",              
#      ylab = "Frequency",          
#      xlim = c(min(size_length$TE_size), max(size_length$TE_size)),
#      breaks = 15)

TE_size_anno <- merge(size_length,TE_anno)

ggplot(size_length)+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "TE Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey() +
  theme_classic()

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsize.png"), plot=ggplot2::last_plot())

L1_size <- ggplot(data=subset(TE_size_anno, ANNO == "LINE/L1"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "L1 Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
Alu_size <- ggplot(data=subset(TE_size_anno, ANNO == "SINE/Alu"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "Alu Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
SVA_size <- ggplot(data=subset(TE_size_anno, ANNO == "Retroposon/SVA"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "SVA Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
ERV_size <- ggplot(data=subset(TE_size_anno, ANNO == "LTR/ERV"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "ERV Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
combined <- L1_size+Alu_size+SVA_size+ERV_size
combined

## Karyoplote

super_freq_7 <- super_freq_7[order(super_freq_7[,1],super_freq_7[,2]),]

super_freq_Gr <- makeGRangesFromDataFrame(df=super_freq_7,
                                          keep.extra.columns=T,
                                          ignore.strand=T,
                                          seqnames.field="CHR",
                                          start.field="POS",
                                          end.field="end",
                                          starts.in.df.are.0based=FALSE)

kpden <- plotKaryotype(genome="hg38", plot.type = 2)
kpPlotDensity(kpden, data=super_freq_Gr,col="goldenrod")
kpAddMainTitle(kpden,"TE density distribution across the chromosomes")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_karyotype.png"), plot=ggplot2::last_plot())

kpli <- plotKaryotype(genome="hg38", plot.type = 2)
kpLines(kpli, data=super_freq_Gr, y=super_freq_Gr$FREQ_ALT,data.panel=1, col="darkgoldenrod",r0=0.25, r1=1)
kpAxis(kpli, ymax=1, r0=0.25, r1=1,numticks = 2, col="#666666", cex=0.5, text.col="black")
kpAddMainTitle(kpli,"Global Allele freqeuncy")

# TE counts per chr -------------------------------------------------------

biTE_count <- super_freq_3_chr_length %>% dplyr::select(c("CHR","POS","TE_size","chr_length"))
biTE_count <- merge(biTE_count,TE_anno)
biTE_count <- biTE_count[order(biTE_count[,1],biTE_count[,2]),]

# no duplicates are present
# dup_indices <- duplicated(biTE_count$POS) &
#   duplicated(biTE_count$POS, fromLast = TRUE) &
#   duplicated(biTE_count$CHR[biTE_count$CHR == "chrX"])
# biTE_count <- biTE_count[!dup_indices, ]

biTE_count_2 <- biTE_count %>%
  dplyr::group_by(CHR, ANNO) %>%
  dplyr::summarize(count = dplyr::n(),
            mean_TE_size = mean(TE_size),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()

biTE_count_2$norm_count <- biTE_count_2$count/biTE_count_2$mean_chr_length

biTE_count_2[biTE_count_2 == "chr1"] <- "chr01"
biTE_count_2[biTE_count_2 == "chr2"] <- "chr02"
biTE_count_2[biTE_count_2 == "chr3"] <- "chr03"
biTE_count_2[biTE_count_2 == "chr4"] <- "chr04"
biTE_count_2[biTE_count_2 == "chr5"] <- "chr05"
biTE_count_2[biTE_count_2 == "chr6"] <- "chr06"
biTE_count_2[biTE_count_2 == "chr7"] <- "chr07"
biTE_count_2[biTE_count_2 == "chr8"] <- "chr08"
biTE_count_2[biTE_count_2 == "chr9"] <- "chr09"

anno_order <- fct_relevel(biTE_count_2$ANNO,"DNA/misc","Retroposon/SVA","LTR/misc","LTR/ERV",
                          "LINE/misc","LINE/L1","SINE/MIR","SINE/Alu")
new_colors <- c("SINE/Alu"="#440154FF", "SINE/MIR"="#453781FF",
                "LINE/L1"="#39558CFF", "LINE/misc"="#238A8DFF",
                "LTR/ERV"="#29AF7FFF","LTR/misc"="#74D055FF",
                "Retroposon/SVA"= "#B8DE29FF","DNA/misc"="#FDE725FF")

ggplot(biTE_count_2, aes(x = (CHR), y=(norm_count), fill = anno_order)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  xlab("Chromosomes")+
  ylab("Normalized TE counts")+
  labs(fill ="TE type")+
  scale_fill_manual(values=c(new_colors))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEcountperchr.png"), plot=ggplot2::last_plot())

contingency_table <- table(biTE_count$CHR)
c_test_result <- chisq.test(contingency_table)
print(c_test_result)
chisq.test(contingency_table)$expected

biTE_count_3 <- aggregate(cbind(count, norm_count) ~ CHR, data = biTE_count_2, sum)

kruskal.test(count ~ CHR, data = biTE_count_3)
kruskal.test(norm_count ~ CHR, data = biTE_count_3)

# Allele freq distribution and pie charts ----------------------------------

## define freq ranges

super_freq_anno$new_class <- lapply(super_freq_anno$FREQ_ALT, function(x) if(x>0.901){
  super_freq_anno$new_class="Fixed"
}else if (0.634 < x && x <= 0.901){
  super_freq_anno$new_class="Major"
}else if (0.367 < x && x <= 0.634){
  super_freq_anno$new_class="Common"
}else if ((0.10) < x && x <= 0.367){
  super_freq_anno$new_class="Rare"
}else 
  super_freq_anno$new_class="Very Rare")
super_freq_anno$new_class <- as.character(unlist(super_freq_anno$new_class))
table(super_freq_anno$new_class)
super_freq_anno$ANNO <- as.character(unlist(super_freq_anno$ANNO))
super_freq_anno_counted <- super_freq_anno %>% dplyr::count(new_class,ANNO)

colnames(super_freq_anno_counted) = c("CLASS","ANNO","NUMBER")
super_freq_anno_counted$CLASS <- as.factor(unlist(super_freq_anno_counted$CLASS))

## freq distribution

anno_order <- fct_relevel(super_freq_anno_counted$ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
new_color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
new_class_order <- fct_relevel(super_freq_anno$new_class,"Fixed","Major","Common","Rare","Very Rare")
super_freq_anno_counted$ANNO <- reorder(super_freq_anno_counted$ANNO, -super_freq_anno_counted$NUMBER)

ggplot(super_freq_anno,mapping=aes(fill = new_class_order))+
  geom_histogram(binwidth = 0.07,aes(x = FREQ_ALT),color="black")+
  labs(
    # title = "Allele frequency distribution",
    x = "TE frequency",
    y = "TE count",
    fill="Frequency ranges")+
  scale_fill_manual(values = color_class) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqspec.png"), plot=ggplot2::last_plot())

p <- ggplot(super_freq_anno_counted, mapping = aes(x=ANNO,
                                                   y=NUMBER,
                                                   fill=anno_order))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#2D718EFF","#29AF7FFF","#56C667FF","#B8DE29FF","#FDE725FF"))+
  # scale_color_manual(values = COLS,breaks=legend_order)+
  labs(
    x = "TE type",
    y = "TE count",
    fill=("TE type")
  ) +
  theme_classic()+
  theme(plot.title = element_text( size = 16, face = "bold", hjust = 0.5))+
  scale_y_continuous(breaks = seq(0,max(super_freq_anno_counted$NUMBER), by = 500))+
  # theme(axis.text.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  facet_grid(~factor(CLASS,levels=c("Very Rare","Rare","Common","Major","Fixed")))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- new_color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqdist.png"), plot=ggplot2::last_plot())

## pie charts 

super_freq_anno_counted_Q4 <- super_freq_anno_counted %>% filter(CLASS == "Fixed")
super_freq_anno_counted_Q3 <- super_freq_anno_counted %>% filter(CLASS == "Major")
super_freq_anno_counted_Q2 <- super_freq_anno_counted %>% filter(CLASS == "Common")
super_freq_anno_counted_Q1 <- super_freq_anno_counted %>% filter(CLASS == "Rare")
super_freq_anno_counted_Q1_rare <- super_freq_anno_counted %>% filter(CLASS == "Very Rare")

A <- ggplot(super_freq_anno_counted_Q4, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Fixed")

B <- ggplot(super_freq_anno_counted_Q3, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Major")

C <- ggplot(super_freq_anno_counted_Q2, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF"))+
  labs(fill="TE type")+
  ggtitle("Common")

D <- ggplot(super_freq_anno_counted_Q1, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Rare")

E <- ggplot(super_freq_anno_counted_Q1_rare, aes(x="", y=NUMBER,
                                                 fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Very Rare")

combined <- E+D+C+B+A & theme(legend.position = "none")
combined + plot_layout(ncol=5, nrow=1, guides = "collect", axis_titles = "collect", tag_level="keep")

# Annovar -----------------------------------------------------------------

## annovar - genomic location count

## out put of cut -f 4 anno_g.hg38_multianno.txt | sort | uniq -c 

annovar_out <- data.frame(
  class = c("downstream","Exons","Intergenes","Introns","ncRNA exons","ncRNA introns","upstream","upstream;downstream","3'UTR","5'UTR"),
  counts = c(39 ,3,3452 ,2390 ,10 ,430 ,34 ,2,41 ,6))
new_class <- annovar_out %>%
  filter(class == 'upstream'| class == "upstream;downstream" | class == "downstream") 
sum(new_class$counts)
upstream_downstream <- data.frame(class = "Upstream/Downstream", counts = sum(new_class$counts))
updated_annovar_out <- rbind(annovar_out, upstream_downstream)
updated_annovar_out <- updated_annovar_out %>%
  filter(class != 'upstream', class != "upstream;downstream" , class != "downstream") 
updated_annovar_out$counts <- as.numeric(updated_annovar_out$counts)
updated_annovar_out$TE <- as.character("TE")
region_order <- fct_relevel(updated_annovar_out$class,"Intergenes","Introns","Exons","5'UTR","3'UTR","Upstream/Downstream","ncRNA introns","ncRNA exons")

ggplot(updated_annovar_out, aes(x=region_order,y=counts,fill=region_order)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = counts), vjust = -0.5, color = "black") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts",
       fill = "Genomic region" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  # + theme(axis.text.x = element_text(angle = 360, hjust = 1))
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## annovar - location per family

annovar_res <- read.table(file = file.path(sys_dir,"PanTE_human/annovar/R_1alt_annovar_out.txt"),
                          header=T)
colnames(annovar_res) <- c("CHR","POS","LOC")
annovar_res$LOC <- gsub("upstream;downstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("upstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("downstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("intergenic", "Intergenes", annovar_res$LOC)
annovar_res$LOC <- gsub("intronic", "Introns", annovar_res$LOC)
annovar_res$LOC <- gsub("exonic", "Exons", annovar_res$LOC)
annovar_res$LOC <- gsub("UTR5", "5'UTR", annovar_res$LOC)
annovar_res$LOC <- gsub("UTR3", "3'UTR", annovar_res$LOC)
annovar_res$LOC <- gsub("ncRNA_intronic", "ncRNA_Introns", annovar_res$LOC)
annovar_res$LOC <- gsub("ncRNA_exonic", "ncRNA_Exons", annovar_res$LOC)
TE_loc <- merge(super_freq_anno,annovar_res)
nrow(TE_loc)
TE_loc_c <- table(TE_loc$LOC, TE_loc$new_class, TE_loc$ANNO) %>% data.frame()
colnames(TE_loc_c) <- c("LOC","Freq","Anno","counts")

region_order <- fct_relevel(TE_loc_c$LOC,"Intergenes","Introns","Exons","5'UTR","3'UTR","Upstream/Downstream","ncRNA_Introns","ncRNA_Exons")
anno_order <- fct_relevel(TE_loc_c$Anno,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
anno_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")
color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
freq_order <- c("Very Rare","Rare","Common","Major","Fixed")
TE_loc_c$Freq <- factor(TE_loc_c$Freq, levels = freq_order)
TE_loc_c$LOC <- reorder(TE_loc_c$LOC, -TE_loc_c$counts)
TE_loc_c$Anno <- reorder(TE_loc_c$Anno, -TE_loc_c$counts)

P <- ggplot(TE_loc_c, aes(x=Anno,y=(log10(counts)),fill=region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts (log10)",
       fill = "Genomic region" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~Freq,nrow=1,ncol=5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- ggplot_gtable(ggplot_build(P))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

p <- ggplot(TE_loc_c, aes(x=Anno,y=(counts),fill=region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts",
       fill = "TE types" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  # scale_fill_manual(values = anno_colors)+
  facet_wrap(~Freq,nrow=1,ncol=5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

## combine freq distribution and annovar

TE_newloc <- TE_loc_c
TE_newloc$LOC <- gsub("5'UTR", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("3'UTR", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("Upstream/Downstream", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("ncRNA_Exons", "Exons", TE_newloc$LOC)
TE_newloc$LOC <- gsub("ncRNA_Introns","Introns", TE_newloc$LOC)

region_order <- fct_relevel(TE_newloc$LOC,"Intergenes","Regulatory","Introns","Exons")
anno_order <- fct_relevel(TE_newloc$Anno,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
anno_order_by_type <- fct_relevel(TE_newloc$Anno, "SINE/Alu", "SINE/MIR", "LINE/L1", "LINE/misc", "LTR/ERV", "LTR/misc", "Retroposon/SVA", "DNA/misc")
anno_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")
color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
freq_order <- c("Very Rare","Rare","Common","Major","Fixed")
TE_newloc$Freq <- factor(TE_newloc$Freq, levels = freq_order)
TE_newloc$LOC <- reorder(TE_newloc$LOC, -TE_newloc$counts)
TE_newloc$Anno <- reorder(TE_newloc$Anno, -TE_newloc$counts)

anno_levels <- c("SINE/Alu", "SINE/MIR", "LINE/L1", "LINE/misc", "LTR/ERV", "LTR/misc", "Retroposon/SVA", "DNA/misc")
color_gradients <- lapply(anno_colors, function(color) {
  colorRampPalette(c(color, "#F4F3EE"))(4)
})
names(color_gradients) <- anno_levels
color_mapping <- unlist(lapply(seq_along(color_gradients), function(i) {
  setNames(color_gradients[[i]], paste0(names(color_gradients)[i], "_", 1:4))
}))
TE_newloc$Color_Group <- paste0(TE_newloc$Anno, "_", as.numeric(region_order))

p <- ggplot(TE_newloc, aes(x=Anno, y= log10(counts+1), fill = Color_Group)) +
  geom_bar(stat = "identity", position = "stack",) +
  labs(title = "",
       x = "TE types",
       y = "log10 (counts+1)",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = color_mapping) +
  facet_wrap(~Freq, nrow = 1, ncol = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqanno.png"), plot=ggplot2::last_plot())

ggplot(TE_newloc, aes(x= Anno, y= log10(counts+1), fill = region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(fill = "Genomic region") +
  theme_classic() +
  scale_fill_manual(values = colorRampPalette(c("#463f3a","#F4F3EE"))(4)) +
  theme_void()

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqannolegand.png"), plot=ggplot2::last_plot())

# p <- ggplot(TE_newloc, aes(x= Anno, y= log10(counts + 1), fill= anno_order))+
#   geom_bar_pattern(stat = "identity",
#                    pattern_color = "white",
#                    pattern_fill = "white",
#                    pattern_density = 0.35,
#                    aes(pattern = LOC, pattern_angle = LOC, pattern_spacing = LOC))+
#   labs(x = "TE types",
#        y = "log10 (counts+1)",
#        fill = "TE types",
#        pattern = "Genomic location")+
#   theme_classic()+
#   scale_pattern_manual(
#     values = c("none","stripe","circle","wave"),
#     guide = guide_legend(override.aes = list(pattern_fill = 'grey')))+
#   scale_fill_manual(values = anno_colors) +
#   facet_wrap(~Freq, nrow = 1, ncol = 5) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# g <- ggplot_gtable(ggplot_build(p))
# striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
# fills <- color_class
# k <- 1
# for (i in striprt) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)

# Recombination rate, TE size, and gene distance ------------------------------------------------------

## define recombination rate date

rec_rate <- read.table(file = file.path(sys_dir,"PanTE_human/recombination_rate/recombAvg.bed"), header = F,
                       sep = "",
                       dec = ".") %>% dplyr::select(c(-4))
colnames(rec_rate) = c("chr","start","end","rate")

## create 10 000 bins for recombination rate

chromo <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
            "chr21", "chr22","chrX")
bin          <- 10000

df2          <- rec_rate %>% 
  mutate(window = start %/% bin) %>% 
  group_by(window,chr) %>%
  summarise(mean = mean(rate)) %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) 

df_out       <- c()
for(i in unique(df2$chr)) {
  df2_tmp      <- df2 %>% filter(chr == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(chr    = i,
           mean = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(chr,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

rec_rate_binned <- df_out %>% arrange(chr, window) 

fc.ymax <- ceiling(max(abs(range(rec_rate_binned$mean))))
fc.ymin <- -fc.ymax
col.over <- "tomato"
col.under <- "gray"
sign.col <- rep(col.over, length(rec_rate_binned))
sign.col[rec_rate_binned$mean>150] <- col.under

ggplot(rec_rate_binned, mapping=aes(x=window,y=mean,col=sign.col))+
  geom_point()+
  facet_wrap(~factor(chr, levels=chromo))

rec_rate_binned_Gr <- makeGRangesFromDataFrame(df=rec_rate_binned,
                                               keep.extra.columns=T,
                                               ignore.strand=T,
                                               seqnames.field="chr",
                                               start.field="start",
                                               end.field="stop",
                                               starts.in.df.are.0based=FALSE)

rec_rate_binned_noNA <- rec_rate_binned %>% drop_na()

ggplot()+
  geom_boxplot(rec_rate_binned_noNA, mapping=aes(x=chr,y=log10(mean)))

## define gene distance and size for super_freq_anno

## size

super_freq_anno$SIZE <- nchar(super_freq_anno$ALT)

## gene distance rGREAT

set.seed(123)
TE_anno_job = submitGreatJob(TE_anno_Gr,species="hg38")
table = getEnrichmentTables(TE_anno_job)
str(table)
plotRegionGeneAssociations(TE_anno_job)
plotRegionGeneAssociations(TE_anno_job,ontology = "GO Biological Process",term_id="GO:0003251")
availableOntologies(TE_anno_job)
res = great(TE_anno_Gr, "MSigDB:H", "txdb:hg38")
plotRegionGeneAssociations(res)
plotVolcano(res)

TE_distTSS <- data.frame(TE_anno_job@association_tables[["all"]][["chr"]],
                           (TE_anno_job@association_tables[["all"]][["start"]]+1),
                           abs(TE_anno_job@association_tables[["all"]][["distTSS"]]),
                           TE_anno_job@association_tables[["all"]][["name"]])
colnames(TE_distTSS) <- c("CHR","POS","distTSS","names")

min_indices <- tapply(seq_along(TE_distTSS$distTSS), TE_distTSS$names, function(x) x[which.min(TE_distTSS$distTSS[x])])
TE_distTSS_filtered <- TE_distTSS[unlist(min_indices), ]

# TE_distTSS <- TE_distTSS[duplicated(TE_distTSS$names),] ######## alawys keep the first one (need to get new gene dis)

TE_info_distTSS <- merge(super_freq_anno,TE_distTSS_filtered)
TE_info_distTSS$ANNO <- as.factor(unlist(TE_info_distTSS$ANNO))
TE_info_distTSS <- TE_info_distTSS[order(TE_info_distTSS[,1],TE_info_distTSS[,2]),]
TE_size_distTSS <- TE_info_distTSS %>%
  dplyr::select(c("CHR","POS","END","FREQ_ALT","distTSS","SIZE","ANNO"))

## create 10 000 bins for TE_size_distTSS

bin          <- 10000

df2          <- TE_size_distTSS %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_anno = list(ANNO),across(c(FREQ_ALT, distTSS,SIZE), mean))
df2[,6:7] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_anno","mean_altfreq","mean_disttss","mean_size","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_anno","mean_altfreq","mean_disttss","mean_size"))

df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_altfreq = NA,
           mean_disttss = NA,
           mean_size = NA,
           sum_anno = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

TE_size_distTSS_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

fc.ymax <- ceiling(max(abs(range(TE_size_distTSS_binned$mean_altfreq))))
fc.ymin <- -fc.ymax
col.over <- "tomato"
col.under <- "gray"
sign.col <- rep(col.over, length(TE_size_distTSS_binned))
sign.col[TE_size_distTSS_binned$mean_altfreq>0.2] <- col.under

ggplot(TE_size_distTSS_binned, mapping=aes(x=window,y=mean_altfreq,col=sign.col))+
  geom_point()+
  facet_wrap(~factor(CHR, levels=chromo))

colnames(TE_size_distTSS_binned)
colnames(rec_rate_binned_noNA) = c("window","CHR","mean_recrate","pos","end")
TE_size_distTSS_rec_rate <- merge(rec_rate_binned_noNA,TE_size_distTSS_binned)
TE_size_distTSS_rec_rate$new_class <- lapply(TE_size_distTSS_rec_rate$mean_altfreq, function(x) if(x>0.901){
  TE_size_distTSS_rec_rate$new_class="Fixed"
}else if (0.634 < x && x <= 0.901){
  TE_size_distTSS_rec_rate$new_class="Major"
}else if (0.367 < x && x <= 0.634){
  TE_size_distTSS_rec_rate$new_class="Common"
}else if ((0.10) < x && x <= 0.367){
  TE_size_distTSS_rec_rate$new_class="Rare"
}else 
  TE_size_distTSS_rec_rate$new_class="Very Rare")
TE_size_distTSS_rec_rate$new_class <- as.character(unlist(TE_size_distTSS_rec_rate$new_class))
table(TE_size_distTSS_rec_rate$new_class)
TE_genomic_windows <- TE_size_distTSS_rec_rate %>% arrange(CHR,window) %>% 
  dplyr::select(c("CHR","start","stop","window","sum_anno","mean_recrate","mean_altfreq","mean_disttss","mean_size","new_class"))
TE_genomic_windows <- TE_genomic_windows[order(TE_genomic_windows[,1],TE_genomic_windows[,2]),]

## TE_genomic_windows is a data frame that has recombination rate, gene distance, size, Global alt allele freq, annotation, and freq ranges information for all TE vars 

## rec rate boxplots, scatter plots, and stats

## filter outliers (SD > 3)

threshold <- 3
z_scores <- scale(TE_genomic_windows$mean_recrate)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_genomic_windows[-outliers, ]

z_scores <- scale(TE_dataframe_clean$mean_disttss)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_dataframe_clean[-outliers, ]

z_scores <- scale(TE_dataframe_clean$mean_size)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_dataframe_clean[-outliers, ]

TE_dataframe <- TE_dataframe_clean %>%
  # filter(mean_recrate < 500, mean_recrate != 0) %>%  # %>% filter(mean_altfreq > 0.1, mean_altfreq < 0.90) # %>% filter(mean_size < 0.98) %>% filter(mean_disttss < 0.98) %>% 
  mutate(freqeuncy = fct_relevel(new_class, "Fixed", "Major","Common","Rare","Very Rare")) %>% 
  filter(mean_recrate > 0) %>%
  filter(mean_altfreq > 0) %>% 
  filter(mean_size >= 100) %>% 
  filter(CHR != "chrX")

## log transformation

TE_dataframe_log <- TE_dataframe
TE_dataframe_log$mean_size <- log10(TE_dataframe_log$mean_size)
TE_dataframe_log$mean_recrate <- log10(TE_dataframe_log$mean_recrate)
TE_dataframe_log$mean_disttss <- log10(TE_dataframe_log$mean_disttss)

color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")

# Boxplots for recombination rate, TE size, and gene distance -----------------------------------------------------------

## all TE boxplot

allTE_size <-  ggplot()+
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  labs(y = "TE size (log10)",
       x = "Allele frequency",
       fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

allTE_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",
       fill = "Frequency ranges"
  )+
scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
freq_order <- fct_relevel(TE_dataframe_log$new_class ,"Very Rare","Rare","Common","Major","Fixed")

allTE_recrate <- ggplot() +
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freq_order,y= mean_recrate,fill=freqeuncy),outlier.shape = "|",outlier.size = 3)+
  labs(y="log10(Recombination rate)",
       x= "TE frequency",
       fill= "Frequency ranges")+
  # geom_jitter(color="black", size=0.2, alpha=0.2)
  scale_fill_manual(values = color_class)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate_boxplot_2.png"), plot=ggplot2::last_plot())

# allTE_size <- ggplot(TE_dataframe_log,aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) +
#   geom_lv(color='black', size=0.75) +
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length (log10 scaled)",
#        x = "Allele frequency"
#        )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("All TE families")+
#   scale_fill_viridis(discrete = TRUE,option="H")+
#   theme_minimal() +
#   theme(
#     panel.grid.major.y = element_line(color = "white", linetype = "dashed"),
#     axis.line = element_line(color = "white"),
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text = element_text(size = 10))
# allTE_disttss <- ggplot(TE_dataframe_log,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene (log10 scaled)",
#      x = "Allele frequency"
#   )+ # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
#   scale_fill_viridis(discrete = TRUE,option="H")+
#   theme_minimal() +
#   theme(
#     panel.grid.major.y = element_line(color = "white", linetype = "dashed"),
#     axis.line = element_line(color = "white"),
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text = element_text(size = 10))

color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
freq_order <- fct_relevel(TE_dataframe_log$new_class ,"Very Rare","Rare","Common","Major","Fixed")

allTE_recrate <-ggplot(TE_dataframe_log,aes(x = freq_order, y=mean_recrate,fill=freqeuncy)) +
  geom_lv(color='black', size=0.75,outlier.shape = "|",outlier.size = 3) +
  geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000",outlier.shape = "|") +
  labs(y="log10(Recombination rate)",
       x= "TE frequency",
       fill= "Frequency ranges")+
  # geom_jitter(color="black", size=0.2, alpha=0.2)
  scale_fill_manual(values = color_class)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate_boxplot.png"), plot=ggplot2::last_plot())

## L1, Alu, SVA boxplot

TE_dataframe_log$new_class <- (unlist(TE_dataframe_log$new_class))
TE_dataframe_log_nolist <- TE_dataframe_log %>%  unnest(sum_anno)

TE_dataframe_log_Alu <- TE_dataframe_log_nolist %>% filter(sum_anno == "SINE/Alu")
TE_dataframe_log_L1 <- TE_dataframe_log_nolist %>% filter(sum_anno == "LINE/L1")
TE_dataframe_log_SVA <- TE_dataframe_log_nolist %>% filter(sum_anno == "Retroposon/SVA")

## SVA

# SVA_size <- ggplot(TE_dataframe_log_SVA, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("SVA")
SVA_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# SVA_disttss <- ggplot(TE_dataframe_log_SVA,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
SVA_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# SVA_recrate <-ggplot(TE_dataframe_log_SVA,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   ) # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
SVA_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+   ggtitle("SVA")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## Alu

# Alu_size <- ggplot(TE_dataframe_log_Alu, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("Alu")
Alu_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Alu_disttss <- ggplot(TE_dataframe_log_Alu,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   ) # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
Alu_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Alu_recrate <-ggplot(TE_dataframe_log_Alu,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
Alu_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  ggtitle("Alu")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## L1

# L1_size <- ggplot(TE_dataframe_log_L1, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("L1")
L1_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# L1_disttss <- ggplot(TE_dataframe_log_L1,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
L1_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# L1_recrate <-ggplot(TE_dataframe_log_L1,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
L1_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+   ggtitle("L1")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

combined <- allTE_recrate + Alu_recrate + L1_recrate + SVA_recrate +
  allTE_disttss + Alu_disttss + L1_disttss + SVA_disttss +
  allTE_size + Alu_size + L1_size + SVA_size &
  theme(legend.position = "bottom")
combined + plot_layout(ncol=4, nrow=3, guides = "collect",
                       axis_titles = "collect")

# Scatter plots for recombination rate, TE size, and gene distance -------------------------------------------------

color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")

a <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_altfreq,x=mean_size,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(y="TE size",
       x= "Allele frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
b <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_altfreq,x=mean_disttss,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(y="Gene distance",
       x= "Allele frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
c <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_altfreq,x=mean_recrate,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(x="Recombination rate",
       y= "TE frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate.png"), plot=ggplot2::last_plot())

combined <- a + b + c & theme(legend.position = "bottom")
combined + plot_layout(ncol=3, nrow=3, guides = "collect",
                       axis_titles = "collect")

combined <- c + b + a +
  allTE_recrate + allTE_disttss + allTE_size +
  Alu_recrate + Alu_disttss + Alu_size + 
  L1_recrate + L1_disttss + L1_size + 
  SVA_recrate + SVA_disttss + SVA_size & theme(legend.position = "bottom")
combined + plot_layout(ncol=3, nrow=5, guides = "collect",
                       axis_titles = "collect")

# Stats and background intervals for recombination rate, TE size, and gene distance ----------------

## correlations

plot(TE_dataframe$mean_size,TE_dataframe$mean_altfreq)
plot(TE_dataframe$mean_disttss,TE_dataframe$mean_altfreq)
plot(TE_dataframe$mean_recrate,TE_dataframe$mean_altfreq)

cor(TE_dataframe$mean_size,TE_dataframe$mean_altfreq)
cor(TE_dataframe$mean_disttss,TE_dataframe$mean_altfreq)
cor(TE_dataframe$mean_recrate,TE_dataframe$mean_altfreq)

pairs(~mean_altfreq+mean_disttss+mean_recrate+mean_size,data=TE_dataframe,
      pch = 21,
      diag.panel=NULL,
      upper.panel = panel.smooth,
      bg = TE_dataframe$sum_anno,
      main = "All TE types"
)

## stats for low/high rec rates regions

TE_dataframe 

add_recrate_category <- function(data, threshold = 5) {
  data <- mutate(data,
                 recrate_category = ifelse(mean_recrate > threshold, "high_recrate", "low_recrate"))
  return(data)
}
TE_dataframe_recratecat <- add_recrate_category(TE_dataframe, threshold = 5)

color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
color_class <- c("#fff3d1","#fed976","#feb24c","#fd8d3c", "#fc4e2a")
freq_order <- fct_relevel(TE_dataframe_recratecat$new_class ,"Very Rare","Rare","Common","Major","Fixed")
recrate_order <- fct_relevel(TE_dataframe_recratecat$recrate_category ,"low_recrate","high_recrate")

ggplot()+
  geom_boxplot(TE_dataframe_recratecat,mapping=aes(x= freq_order,y=log10(mean_recrate),fill=recrate_order),outlier.shape = "|",outlier.size = 3)+
  labs(y="log10(Recombination rate)",
       x= "TE frequency",
       fill= "Recombination rate")+
  scale_fill_manual(values = c("#4B878BFF","#D01C1FFF"))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recrate_cat.png"), plot=ggplot2::last_plot())

TE_recrate <- table(TE_dataframe_recratecat$recrate_category,TE_dataframe_recratecat$new_class) %>% as.data.frame()
colnames(TE_recrate) <- c("recrate_category","new_class","count")
freq_order <- fct_relevel(TE_recrate$new_class ,"Very Rare","Rare","Common","Major","Fixed")
recrate_order <- fct_relevel(TE_recrate$recrate_category ,"low_recrate","high_recrate")
ggplot(TE_recrate,mapping=aes(x=freq_order,y=count,fill=recrate_category))+
  geom_bar(position="fill", stat="identity")+
  labs(y="Count",
       x= "TE frequency",
       fill= "Recombination rate")+
  scale_fill_manual(values = c("#D01C1FFF","#4B878BFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recrate_cat_2.png"), plot=ggplot2::last_plot())

low_df <- TE_recrate %>% filter(recrate_category == "low_recrate")
high_df <- TE_recrate %>% filter(recrate_category == "high_recrate")
categories <- unique(TE_recrate$new_class)
results <- lapply(categories, function(category) {
  low_count <- low_df$count[low_df$new_class == category]
  high_count <- high_df$count[high_df$new_class == category]
  chi_squared <- chisq.test(c(low_count, high_count))
  list(
    category = category,
    chi_squared = chi_squared$statistic,
    df = chi_squared$parameter,
    p_value = chi_squared$p.value
  )
})
for (result in results) {
  cat("Category:", result$category, "\n")
  cat("Chi-squared value:", result$chi_squared, "\n")
  cat("Degrees of freedom:", result$df, "\n")
  cat("P-value:", result$p_value, "\n\n")}

## percentage of counts in low rec region:

(3166/(3166+547))*100

## background intervals

background_hg38 <- read.table(file = file.path(sys_dir,"PanTE_human/recombination_rate/background_3713_hg38.bed"),header=F)
colnames(background_hg38) <- c("CHR","POS","END","id","size","strand")
background_hg38 <- background_hg38[order(background_hg38[,1],background_hg38[,2]),]
background_hg38_Gr <- makeGRangesFromDataFrame(df=background_hg38,
                                               keep.extra.columns=T,
                                               ignore.strand=F,
                                               seqnames.field="CHR",
                                               start.field="POS",
                                               end.field="END",
                                               starts.in.df.are.0based=FALSE)
set.seed(123)
background_hg38_job = submitGreatJob(background_hg38_Gr,species="hg38")
plotRegionGeneAssociations(background_hg38_job)
background_hg38_distTSS <- data.frame(background_hg38_job@association_tables[["all"]][["chr"]],
                                      (background_hg38_job@association_tables[["all"]][["start"]]+1),
                                      abs(background_hg38_job@association_tables[["all"]][["distTSS"]]),
                                      background_hg38_job@association_tables[["all"]][["name"]])
colnames(background_hg38_distTSS) <- c("CHR","POS","distTSS","names")
min_indices <- tapply(seq_along(background_hg38_distTSS$distTSS), background_hg38_distTSS$names, function(x) x[which.min(background_hg38_distTSS$distTSS[x])])
background_hg38_distTSS_filtered <- background_hg38_distTSS[unlist(min_indices), ]

bin          <- 10000

df2          <- background_hg38_distTSS_filtered %>% 
  mutate(window = POS %/% bin) %>% 
  group_by(window,CHR) %>%
  summarise(across(c(distTSS), mean)) %>% 
  mutate(start = (window*bin)+1, stop=(window+1)*bin)
colnames(df2) <- c("window","CHR","mean_distTSS","start","stop")
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_distTSS = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
background_hg38_distTSS_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

threshold <- 3
z_scores <- scale(rec_rate_binned_noNA$mean_recrate)
outliers <- which(abs(z_scores) > threshold)
rec_rate_binned_noNA_filterd <- rec_rate_binned_noNA[-outliers, ]
background_rec <- merge(rec_rate_binned_noNA,background_hg38_distTSS_binned)
background_rec <- merge(rec_rate_binned_noNA_filterd,background_hg38_distTSS_binned) %>% drop_na() %>% 
  dplyr::select(c("CHR","start","stop","window","mean_recrate","mean_distTSS"))
colnames(background_rec) <- (c("CHR","start","stop","window","mean_recrate_background","mean_distTSS"))

table(TE_dataframe$new_class)
# background_rec_less100 <- background_rec %>% filter(mean_recrate_background < 100)
# TE_dataframe_less100 <- TE_dataframe %>% filter(mean_recrate < 100)

plot(data = background_rec, mean_recrate_background ~ window, col="lightgray")
plot(data = TE_dataframe, mean_recrate ~ window, col="#EE353E")

ggplot()+
  geom_point(data=subset(background_rec,mean_recrate_background < 1000),mapping=aes(x=window,y=mean_recrate_background),shape=25,col="lightgray")+
  geom_point(data=subset(TE_dataframe,mean_recrate < 1000),mapping=aes(x=window,y=mean_recrate),shape=25,col="#EE353E",fill="#EE353E")+
  theme_classic()+
  labs(
    x="Genomic window",
    y="Recombination rate")

ggplot()+
  geom_point(data=subset(background_rec,mean_recrate_background < 50),mapping=aes(x=window,y=mean_recrate_background),shape=25,col="lightgray")+
  geom_point(data=subset(TE_dataframe,mean_recrate < 50),mapping=aes(x=window,y=mean_recrate),shape=25,col="#EE353E",fill="#EE353E")+
  theme_classic()+
  labs(
    x="Genomic window",
    y="Recombination rate")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recratebackground.png"), plot=ggplot2::last_plot())

background_rec <- background_rec %>% mutate(Source = "Background")
TE_dataframe <- TE_dataframe %>% mutate(Source = "TE")

ggplot()+
  # geom_jitter(data = background_rec, aes(x = CHR, y = mean_recrate_background,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(background_rec,mapping=aes(x=CHR,y=(mean_recrate_background),group=Source),fill=alpha('black',0.5),col="black",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))
  
ggplot()+
  # geom_jitter(data = TE_dataframe, aes(x = CHR, y = mean_recrate,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(TE_dataframe,mapping=aes(x=CHR,y=(mean_recrate),group=Source),fill=alpha('#EE353E',0.5),col="#EE353E",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))

wilcox_test <- wilcox.test(TE_dataframe$mean_recrate,background_rec$mean_recrate_background)
print(wilcox_test)
length(TE_dataframe$mean_recrate)
length(background_rec$mean_recrate_background)

# TE divergence in Global, AFR, and outofAFR -----------------------------------------------------------

## define kimura distance from RepeatMasker

panic_KD <- read.csv(file.path(sys_dir,"PanTE_human/evo_age/Kdistance_filtered.txt"),
                     header = F,
                     sep="",
                     dec = ".")
extract_last_part <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) > 1 && grepl("Kimura", x)) {
    return(as.numeric(parts[length(parts)]))
  } else {
    return(x)
  }
}
panic_KD$chr_div <- sapply(panic_KD$V1, extract_last_part)
chromosome_indices <- grep("chromosome", panic_KD$chr_div)
for (i in 1:(length(chromosome_indices) - 1)) {
  start_index <- chromosome_indices[i] + 1
  end_index <- chromosome_indices[i + 1] - 1
  if (end_index > start_index) {
    values <- as.numeric(panic_KD$chr_div[start_index:end_index])
    mean_value <- mean(values, na.rm = TRUE)
    panic_KD$chr_div[start_index:end_index] <- mean_value
  }
}
panic_KD <- panic_KD %>% dplyr::select(c(-1))
for (i in 1:(nrow(panic_KD) - 1)) {
  if (panic_KD$chr_div[i] == panic_KD$chr_div[i + 1]) {
    panic_KD$chr_div[i] <- NA
  }
}
panic_KD <- panic_KD %>% drop_na()
panic_KD <- panic_KD[-10361, , drop = FALSE]
panic_KD$div <- NA
panic_KD$chr <- NA
panic_KD$div[seq(2, nrow(panic_KD), by = 2)] <- panic_KD$chr_div[seq(2, nrow(panic_KD), by = 2)]
panic_KD$chr[seq(1, nrow(panic_KD), by = 2)] <- panic_KD$chr_div[seq(1, nrow(panic_KD), by = 2)]
panic_KD_div <- data.frame(panic_KD$div) %>% drop_na()
panic_KD_chr <- data.frame(panic_KD$chr) %>% drop_na()
length(panic_KD_div$panic_KD.div)
length(panic_KD_chr$panic_KD.chr)
panic_KD_2 <- data.frame(cbind(panic_KD_div$panic_KD.div,panic_KD_chr$panic_KD.chr))
colnames(panic_KD_2) = c("div","chr")
panic_KD_2$chr <- gsub("_", "\t", panic_KD_2$chr)
panic_KD_2[,3:5] <- stringr::str_split_fixed(panic_KD_2$chr, "\t", 3) 
panic_KD_2 <- panic_KD_2 %>% dplyr::select(c(-"chr"))
colnames(panic_KD_2) = c("div","CHR","POS","svid")
panic_KD_2$svid <- gsub("\t", ">", panic_KD_2$svid)
panic_KD_2$div <- as.numeric(unlist(panic_KD_2$div))
panic_KD_2$POS <- as.numeric(unlist(panic_KD_2$POS))
panic_KD_2$CHR <- gsub("chromosome","chr",panic_KD_2$CHR)
panic_KD_2 <- panic_KD_2 %>% dplyr::select(c("CHR","POS",-"svid","div"))
panic_KD_2 <- panic_KD_2[order(panic_KD_2[,1],panic_KD_2[,2]),]

TE_size_distTSS_div <- merge(TE_size_distTSS,panic_KD_2)
TE_size_distTSS_div$ANNO <- as.character(unlist(TE_size_distTSS_div$ANNO))
TE_size_distTSS_div$div <- as.numeric(TE_size_distTSS_div$div)

## include AFR and outofAFR info

AFR_freq_edited <- AFR_freq_3 %>% dplyr::select(c("CHR","POS","FREQ_ALT"))
colnames(AFR_freq_edited) = c("CHR","POS","FREQ_ALT.a")
outofAFR_freq_edited <- outofAFR_freq_3 %>% dplyr::select(c("CHR","POS","FREQ_ALT.o"))
TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_size_distTSS_div_AO <- merge(TE_africa_outofafrica,TE_size_distTSS_div)

TE_size_distTSS_div_AO$ANNO <- as.character(unlist(TE_size_distTSS_div_AO$ANNO))
TE_size_distTSS_div_AO$div <- as.numeric(TE_size_distTSS_div_AO$div)

## neutral evo age calculation for TE_size_distTSS_div_AO

TE_div_Ao_neu <- TE_size_distTSS_div_AO %>%
  filter(FREQ_ALT < 1)

meo <- (0.016 + 0.025 + 0.016)/3
TE_div_Ao_neu$div_age <- ((TE_div_Ao_neu$div)/(2*(meo)))
Ne <- 12500
x <- TE_div_Ao_neu$FREQ_ALT
TE_div_Ao_neu$exp_age_ancestral <- (-4*Ne*(x/(1-x))*log(x))
Ne <- 24000
x <- TE_div_Ao_neu$FREQ_ALT.a
TE_div_Ao_neu$exp_age_africa <- (-4*Ne*(x/(1-x))*log(x))
Ne <- 7700
x <- TE_div_Ao_neu$FREQ_ALT.o
TE_div_Ao_neu$exp_age_outofafrica <- (-4*Ne*(x/(1-x))*log(x))

get_custom_color <- function(value) {
  if (value > +3) {
    return("Old")
  } else if (value < -3) {
    return("Young")
  } else {
    return("Neutral")  
  }
}

TE_neu_L1AluSVA <- TE_div_Ao_neu %>% 
  filter(ANNO == "SINE/Alu" | ANNO == "LINE/L1" | ANNO == "Retroposon/SVA") # %>% filter(div > 0)
TE_neu_L1AluSVA$diff_age <- scale(scale(TE_neu_L1AluSVA$div_age))-(scale(TE_neu_L1AluSVA$exp_age_ancestral))
percentile_99 <- quantile(TE_neu_L1AluSVA$diff_age, 0.99)
percentile_01 <- quantile(TE_neu_L1AluSVA$diff_age, 0.01)
TE_neu_L1AluSVA$custom_color <- sapply(TE_neu_L1AluSVA$diff_age, get_custom_color)

ggplot()+
  geom_point(data= TE_neu_L1AluSVA,mapping=aes(y=diff_age,x=POS,col=custom_color),size = 1)+
  labs(
    y="Scaled age estimate difference",
    col = "TE age") +
  scale_color_manual(values = c("Old" = "orange", "Young" = "steelblue", "Neutral" = "gray")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_evoagediff.png"), plot=ggplot2::last_plot())

## TE divergence/freq plots - similarity based AFR and outofAFR

anno_order <- c("SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
new_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")

TE_div_Ao_neu$anno_order <- ordered(TE_div_Ao_neu$ANNO,levels = anno_order)
TE_div_Ao_neu$diff_age <- scale(scale(TE_div_Ao_neu$div_age))-(scale(TE_div_Ao_neu$exp_age_ancestral))
TE_div_Ao_neu$custom_color <- sapply(TE_div_Ao_neu$diff_age, get_custom_color)

plot_all <- ggplot()+
  geom_point(TE_div_Ao_neu, mapping=aes(x=scale(div_age),y=FREQ_ALT,col=custom_color))+
  geom_line(TE_div_Ao_neu, mapping=aes(x=scale(exp_age_ancestral),y=FREQ_ALT,col="black"))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
plot_all <- plot_all + scale_colour_manual(values=c("black","grey","orange","steelblue"))
plot_all <- plot_all + labs(x="Divergence age (scaled)",y="Global TE frequency",col="TE type")
plot_all <- plot_all + theme_linedraw()
plot_all +  scale_size(trans ="reverse")
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv.png"), plot=ggplot2::last_plot())

## include TE freq ranges to TE_div_Ao_neu

TE_div_Ao_neu_freq <- TE_div_Ao_neu

TE_div_Ao_neu_freq$new_class <- lapply(TE_div_Ao_neu_freq$FREQ_ALT, function(x) if(x>0.901){
  TE_div_Ao_neu_freq$new_class="Fixed"
}else if (0.634 < x && x <= 0.901){
  TE_div_Ao_neu_freq$new_class="Major"
}else if (0.367 < x && x <= 0.634){
  TE_div_Ao_neu_freq$new_class="Common"
}else if ((0.10) < x && x <= 0.367){
  TE_div_Ao_neu_freq$new_class="Rare"
}else 
  TE_div_Ao_neu_freq$new_class="Very Rare")
TE_div_Ao_neu_freq$new_class <- as.character(unlist(TE_div_Ao_neu_freq$new_class))
table(TE_div_Ao_neu_freq$new_class)
color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
class_order <- fct_relevel(TE_div_Ao_neu_freq$new_class,"Very Rare","Rare","Common","Major","Fixed")

ggplot()+
  geom_point(TE_div_Ao_neu_freq, mapping=aes(x=scale(div_age),y=FREQ_ALT,col=class_order))+
  geom_line(TE_div_Ao_neu_freq, mapping=aes(x=scale(exp_age_ancestral),y=FREQ_ALT,col="black"))+
  scale_colour_manual(values=c(color_class,"black"))+
  labs(x="Divergence age (scaled)",y="Global TE frequency",col="TE type")+
  theme_linedraw()+
  scale_size(trans ="reverse")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv_freqrange.png"), plot=ggplot2::last_plot())

## include candi_all_2 from the following selection scans sections to label candidate TE loci under positive selection
## include TE_Global_shared from the following venn diagram sections to extract shared TE loci

TE_candi_all <- TE_div_Ao_neu

TE_candi_all$candilabel[TE_candi_all$POS %in% candi_all_2$POS] <- "candidate"
TE_candi_all$candilabel[is.na(TE_candi_all$candilabel)] <- "not_candidate"

TE_candi_all$globalsharedlabel[TE_candi_all$POS %in% TE_Global_shared$POS] <- "globalshared"
TE_candi_all$globalsharedlabel[is.na(TE_candi_all$globalsharedlabel)] <- "not_globalshared"

plot_africa <- ggplot()+
  geom_point(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(div_age),y=FREQ_ALT.a,col=candilabel)) +
  geom_line(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(exp_age_ancestral),y=FREQ_ALT,col="black"))+
  # scale_colour_manual(values=c(new_colors,"black"))+
  scale_color_manual(values=c("black","#EE353E","#F5F5F5"))+
  labs(x="Divergence age (scaled)",y="Africa TE frequency",col="TE type")+
  theme_linedraw()+
  scale_size(trans ="reverse")+
  theme(panel.border = element_rect(color = "#00539CFF", size = 1))

plot_outofafrica <- ggplot()+
  geom_point(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(div_age),y=FREQ_ALT.o,col=candilabel)) +
  geom_line(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(exp_age_ancestral),y=FREQ_ALT,col="black"))+
  # scale_colour_manual(values=c(new_colors,"black"))+
  scale_color_manual(values=c("black","#EE353E","#F5F5F5"))+
  labs(x="Divergence age (scaled)",y="Out-of-Africa TE frequency",col="TE type")+
  theme_linedraw()+
  scale_size(trans ="reverse")+
  theme(panel.border = element_rect(color = "#EEA47FFF", size = 1))

combined <- plot_africa + plot_outofafrica &
  theme(legend.position = "bottom")
combined + plot_layout(ncol=2, nrow=1, guides = "collect",
                       axis_titles = "collect")

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv_perpop.png"), plot=ggplot2::last_plot())

## TE divergence/fst plots

Fst <- read.csv(file.path(sys_dir,"PanTE_human/selection_scans/fst/out.weir.fst"), header = T,sep = "")   
colnames(Fst) <- c("CHR","POS","Fst")
TE_div_Fst <- merge(TE_div_Ao_neu,Fst)

TE_div_Fst_filtered <- TE_div_Fst %>% filter(Fst > 0)
Fst_percentile_99 <- quantile(TE_div_Fst_filtered$Fst, 0.99)

get_custom_fst <- function(value_1,value_2) {
  if (value_1 > Fst_percentile_99 && value_2 == "Young")
  {return("Sig")}
  else {
    return("not_sig")
  }}
TE_div_Fst_filtered$custom_fst <- apply(TE_div_Fst_filtered, 1, function(row) get_custom_fst(row["Fst"],row["custom_color"]))

plot_fst <- ggplot()+
  geom_point(data=(TE_div_Fst_filtered), mapping=aes(x=scale(div_age),y=Fst,col=custom_fst))
  # geom_line(data=subset(TE_div_Fst_filtered, FREQ_ALT < 0.5), mapping=aes(x=scale(exp_age_ancestral),y=FREQ_ALT,col="black"))
plot_fst <- plot_fst + scale_colour_manual(values=c("grey","orange","steelblue"))
plot_fst <- plot_fst + scale_colour_manual(values=c("grey","#EE353E"))
plot_fst <- plot_fst + labs(x="Divergence age (scaled)",y="Fst",col="TE type")
plot_fst <- plot_fst + theme_linedraw()
plot_fst +  scale_size(trans ="reverse") 
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_fstdivplot.png"), plot=ggplot2::last_plot())

## TE_div_Fst_filtered is a data frame that has gene distance, size, global, AFR, and outofAFR alt allele freq, annotation, divergence, neutral evo calculation, and Fst information for all TE vars 

candidate_TE_div_Fst <- TE_div_Fst_filtered %>%
  filter(custom_color == "Young") %>%
  filter(Fst > Fst_percentile_99)

# TE_FULL_INFO ------------------------------------------------------------

## create 10 000 bins for TE_div_Fst_filtered

TE_div_Fst_filtered$POS <- as.numeric(TE_div_Fst_filtered$POS)

bin          <- 10000

df2          <- TE_div_Fst_filtered %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_anno = list(ANNO),across(c(FREQ_ALT,FREQ_ALT.a,FREQ_ALT.o,distTSS,SIZE,div,div_age,exp_age_ancestral,diff_age,Fst), mean))
df2[,13:14] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_anno","mean_altfreq","FREQ_ALT.a","FREQ_ALT.o","mean_disttss","mean_size","div","div_age","exp_age_ancestral","diff_age","Fst","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_anno","mean_altfreq","FREQ_ALT.a","FREQ_ALT.o","mean_disttss","mean_size","div","div_age","exp_age_ancestral","diff_age","Fst"))
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_altfreq = NA,
           FREQ_ALT.a = NA,
           FREQ_ALT.o = NA,
           mean_disttss = NA,
           mean_size = NA,
           div = NA,
           div_age = NA,
           exp_age_ancestral = NA,
           diff_age = NA,
           Fst = NA,
           sum_anno = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

TE_div_Fst_binned <- df_out %>% arrange(CHR, window) %>% drop_na()
colnames(TE_div_Fst_binned)
colnames(rec_rate_binned_noNA) = c("window","CHR","mean_recrate","pos","end")
TE_div_Fst_rec_rate <- merge(rec_rate_binned_noNA,TE_div_Fst_binned)
TE_div_Fst_rec_rate$new_class <- lapply(TE_div_Fst_rec_rate$mean_altfreq, function(x) if(x>0.901){
  TE_div_Fst_rec_rate$new_class="Fixed"
}else if (0.634 < x && x <= 0.901){
  TE_div_Fst_rec_rate$new_class="Major"
}else if (0.367 < x && x <= 0.634){
  TE_div_Fst_rec_rate$new_class="Common"
}else if ((0.10) < x && x <= 0.367){
  TE_div_Fst_rec_rate$new_class="Rare"
}else 
  TE_div_Fst_rec_rate$new_class="Very Rare")
TE_div_Fst_rec_rate$new_class <- as.character(unlist(TE_div_Fst_rec_rate$new_class))
table(TE_div_Fst_rec_rate$new_class)
TE_FULL_INFO <- TE_div_Fst_rec_rate %>% arrange(CHR,window) %>% 
  dplyr::select(c("CHR","start","stop","window","mean_recrate","mean_altfreq","FREQ_ALT.a","FREQ_ALT.o","mean_disttss","mean_size",
                  "div","div_age","exp_age_ancestral","diff_age","Fst"))
TE_FULL_INFO <- TE_FULL_INFO[order(TE_FULL_INFO[,1],TE_FULL_INFO[,2]),]

## TE_FULL_INFO is a data frame that has ALL information for all TE vars ## check duplicated entries


# Figure 2 ----------------------------------------------------------------



# TE counts per individual ---------------------------------------------------------------

## define population info, major, singletons, and shared TEs

## major TE

population_info <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/hprc_year1_sample_metadata.txt"), header = F,
                              sep = "",
                              dec = ".")

major_TE_filtered_counted <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/major_TE_stat"), header = F,
                                        sep = "",
                                        dec = ".")
major_filtered_counted_1 <- subset(major_TE_filtered_counted, select=(c(3,9)))
colnames(major_filtered_counted_1) <- c("V1","V2")
maj_pop <- right_join(major_filtered_counted_1, population_info, by='V1') %>%
  arrange(desc(V3),(V2.x)) %>%
  mutate(V1 = factor(V1, levels = V1))
ggplot(maj_pop,aes(x = V1, y = V2.x, fill = V3)) +
  geom_col()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  xlab("ind")+
  ylab("variant counts")+
  ggtitle("major_TE_filtered_counted")+
  labs(fill = "superpop")

## singleton TE

singletons_TE_filtered_counted <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/R_singletons_TE_filtered.vcf.singletons"), header = F,
                                             sep = "",
                                             dec = ".")
col_order <- c("V2","V1")
singletons_TE_filtered_counted <- singletons_TE_filtered_counted[,col_order]
colnames(singletons_TE_filtered_counted) <- c("V1","V2")
sin_pop <- right_join(singletons_TE_filtered_counted, population_info, by='V1') %>%
  arrange(desc(V3),(V2.x)) %>%
  mutate(V1 = factor(V1, levels = V1))
ggplot(sin_pop,aes(x = V1, y = V2.x, fill = V3)) +
  geom_col()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  xlab("ind")+
  ylab("variant counts")+
  ggtitle("singletons_TE_filtered_counted")+
  labs(fill = "superpop")
ggplot(sin_pop,aes(x = reorder(V2.y,V2.x), y = V2.x, fill = V3)) +
  geom_boxplot()+
  theme_minimal()+
  geom_jitter(color="black", size=0.5, alpha=0.9)+
  xlab("subpop")+
  ylab("variant counts")+
  ggtitle("singletons_TE_filtered_counted")+
  labs(fill = "superpop")+
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest1"))

## major + singletons + shared

colnames(sin_pop)=c("ind","singletons","subpop","superpop")
colnames(maj_pop)=c("ind","major","subpop","superpop")
maj_sin <- merge(maj_pop,sin_pop)
maj_sin_shar <- data.frame(append(maj_sin, c(V5=105), after=5))
## counts of shared var = 105. The number might change depending on the nonref TE var calculated in source

colnames(maj_sin_shar)=c("ind","subpop","superpop","Major allele frequency > 0.5","Singletons","Shared in all but GRCh38")
maj_sin_shar <- cbind(maj_sin_shar, Total = rowSums(maj_sin_shar[,c(4,5,6)]))
maj_sin_shar_tidy <- pivot_longer(maj_sin_shar, cols = (c(4,5,6)), values_drop_na = TRUE)
colnames(maj_sin_shar_tidy)=c("ind","subpop","superpop","Total","filter_type","filter_type_count")
maj_sin_shar_tidy <- maj_sin_shar_tidy %>% drop_na() %>% as.data.frame() 
maj_sin_shar_tidy <- maj_sin_shar_tidy %>% filter(maj_sin_shar_tidy$ind != "NA19240")
maj_sin_shar_tidy$filter_type <- as.factor(maj_sin_shar_tidy$filter_type)
maj_sin_shar_tidy$filter_type <- fct_relevel(maj_sin_shar_tidy$filter_type, "Major allele frequency > 0.5","Singletons","Shared in all but GRCh38")

ggplot(maj_sin_shar_tidy,aes(x = reorder(ind,filter_type_count), y = filter_type_count, fill = filter_type)) +
  geom_bar(position="stack", stat="identity",width=0.7)+
  # theme_set(theme_minimal(base_family = "Georgia"))+
  scale_y_continuous(breaks = seq(0, 465, by = 50)) +
  labs(fill = "TEs",x="Individuals",y="TE counts")+
  scale_fill_manual(values =  c("#5F9E9E","#272772","#B22224"))+
  theme_classic()+
  theme()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperind.png"), plot=ggplot2::last_plot())

# TE counts per population and Venn diagram ------------------------------------------------

## TE counts per population

AFR <- c("HG01891","HG02257","HG02486","HG02559","HG02572","HG03516","HG02622","HG02630","HG02717","HG02886","HG03453","HG03540","HG03579","HG02109","HG02145","HG02723","HG02818","HG03486","NA18906","NA20129","NA21309","HG02055","HG03098")
outofAFR <- c("HG01123","HG01258","HG01358","HG01361","HG00735","HG00741","HG01071","HG01106","HG01175","HG01928","HG01952","HG01978","HG02148","HG00733","HG01109","HG01243","HG00438","HG00621","HG00673","HG02080","HG03492")

TE_indv <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/TE_counts_1"),header = F,sep = "") %>%
  dplyr::select(c(3,9)) %>%
  filter(V3 != "CHM13")
colnames(TE_indv) <- c("indv","counts")

pop <- rep(NA, length(TE_indv$indv))
for (value in AFR) {pop[grep(value, TE_indv$indv)] <- "Africa"}
for (value in outofAFR) {pop[grep(value, TE_indv$indv)] <- "outofAfrica"}
TE_indv <- as.tibble(data.frame(TE_indv, pop)) %>% as.data.frame()

ggplot(TE_indv,aes(x = reorder(indv,counts), y = counts,fill=pop)) +
  geom_bar(position="stack", stat="identity",width=0.7)+
  scale_y_continuous(breaks = seq(0,700, by = 50)) +
  labs(fill="Population",x="Individuals",y="TE counts")+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
  )+
  scale_fill_manual(values =  c("#00539CFF","#EEA47FFF"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperindperpop.png"), plot=ggplot2::last_plot())

TE_indv_AFR <- TE_indv %>% filter(pop=="Africa")
TE_indv_outofAFR <- TE_indv %>% filter(pop=="outofAfrica")
summary(TE_indv_AFR$counts)
sd(TE_indv_AFR$counts)
summary(TE_indv_outofAFR$counts)
sd(TE_indv_outofAFR$counts)

shapiro.test(TE_indv$counts[TE_indv$pop == "Africa"])
shapiro.test(TE_indv$counts[TE_indv$pop == "outofAfrica"])

wilcox_test_result <- wilcox.test(counts ~ pop, data = TE_indv)
print(wilcox_test_result)

p <- ggplot(TE_indv,aes(x = pop, y = counts, fill=pop)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),color="black",alpha=0.3)+
  # geom_text(aes(x = 1.5, y = max(TE_indv$counts)+10, label = significance_level), size = 6)+
  scale_y_continuous(breaks = seq(0,700, by = 50)) +
  labs(fill="Population",x="Population",y="TE counts")+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )+
  # annotate("text", x = 1.5, y = 650, label = paste("Wilcoxon, p-value=", round(wilcox_test_result$p.value,8)), size = 4)+
  scale_fill_manual(values =  c("#00539CFF","#EEA47FFF"))+
  theme_classic()
p + stat_compare_means(method = "wilcox", label.y = 680)

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperpopboxplot.png"), plot=ggplot2::last_plot())

## Venn diagram 

# TE_africa_outofafrica_no_dup <- TE_africa_outofafrica %>% filter(!duplicated(.))

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_africa_outofafrica$id <- seq(1,6407)

TE_OO <-TE_africa_outofafrica %>% filter(FREQ_ALT.o == 0)
TE_AA <-TE_africa_outofafrica %>% filter(FREQ_ALT.a == 0)
TE_O <-TE_africa_outofafrica %>% filter(FREQ_ALT.o > 0) 
TE_A <- TE_africa_outofafrica %>% filter(FREQ_ALT.a > 0) 
TE_Global <- TE_africa_outofafrica
TE_Global_shared <- TE_africa_outofafrica %>%
  filter(FREQ_ALT.a != 0) %>%
  filter(FREQ_ALT.o != 0) 
TE_Global_zero <- TE_africa_outofafrica %>%
  filter(FREQ_ALT.a == 0 & FREQ_ALT.o == 0) 

a <- length(TE_OO$id)
b <- length(TE_AA$id)
c <- length(TE_O$id)
d <- length(TE_A$id)
e <- length(TE_Global$id)
f <- length(TE_Global_shared$id)
g <- length(TE_Global_zero$id)

venn.plot <- venn.diagram(x= list(TE_A$id, TE_O$id),
                          filename = "fig2_vendia.png",
                          category.names = c("Africa","Out-of-Africa"),
                          output=TRUE,
                          col=c("#00539CFF", '#EEA47FFF'),
                          fill = c(alpha("#00539CFF",0.3), alpha('#EEA47FFF',0.3)),
                          fontface = "bold",
                          cat.fontface = "bold",
                          cat.cex = 0.4,
                          cat.pos = c(-27, 27),
                          cat.default.pos = "outer",
                          cex = 0.6,
                          # cat.pos = c(0, 0),
                          cat.dist = c(0.05, 0.05),
                          imagetype="png" ,
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw")

# Allele frequency per population -----------------------------------------

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
# TE_africa_outofafrica_no_dup <- do.call(rbind, lapply(split(TE_africa_outofafrica, TE_africa_outofafrica$CHR), function(group) group[!duplicated(group), ]))
# unique(length(TE_africa_outofafrica_no_dup$POS))
TE_africa_outofafrica$id <- seq(1,6407)
TE_africa_outofafrica$pop <- c("population")
TE_africa <- TE_africa_outofafrica %>% dplyr::select(c(1,2,3,5)) %>% filter(FREQ_ALT.a > 0)
TE_africa$pop <- c("Africa")
TE_outofafrica <- TE_africa_outofafrica %>% dplyr::select(c(1,2,4,5)) %>% filter(FREQ_ALT.o > 0)
TE_outofafrica$pop <- c("Out of Africa")

ggplot()+
  geom_jitter(data = TE_outofafrica, aes(x = pop, y = FREQ_ALT.o), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_jitter(data = TE_africa, aes(x = pop, y = FREQ_ALT.a), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(TE_outofafrica,mapping=aes(x=pop,y=(FREQ_ALT.o)),fill=alpha('#EEA47FFF',0.5),col="#EEA47FFF",outlier.shape = "|")+
  geom_boxplot(TE_africa,mapping=aes(x=pop,y=(FREQ_ALT.a)),fill=alpha('#00539CFF',0.5),col="#00539CFF",outlier.shape = "|")+
  labs(x="",y="TE frequency")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold",color="black"))

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_freqperpopbox.png"), plot=ggplot2::last_plot())

# TE family distribution per population - pie charts ----------------------

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_africa_outofafrica$id <- seq(1,6407)
TE_africa_outofafrica$pop <- c("population")
TE_africa_outofafrica_anno <- merge(TE_africa_outofafrica,TE_anno) %>%
  dplyr::select(c("CHR","POS","FREQ_ALT.a","FREQ_ALT.o","id","ANNO"))
TE_africa_anno <- TE_africa_outofafrica_anno %>% dplyr::select(c(1,2,3,5,6)) %>% filter(FREQ_ALT.a > 0)
TE_africa_anno$pop <- c("Africa")
TE_outofafrica_anno <- TE_africa_outofafrica_anno %>% dplyr::select(c(1,2,4,5,6)) %>% filter(FREQ_ALT.o > 0)
TE_outofafrica_anno$pop <- c("Out of Africa")

TE_africa_anno_counted <- TE_africa_anno %>% dplyr::count(pop,ANNO)
colnames(TE_africa_anno_counted) = c("pop","ANNO","NUMBER")
TE_outofafrica_anno_counted <- TE_outofafrica_anno %>% dplyr::count(pop,ANNO)
colnames(TE_outofafrica_anno_counted) = c("pop","ANNO","NUMBER")

AFRpie <- ggplot(TE_africa_anno_counted, aes(x="", y=NUMBER,
                                             fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("AFR")

outofAFRpie <- ggplot(TE_outofafrica_anno_counted, aes(x="", y=NUMBER,
                                                       fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("outofAFR")

combined <- AFRpie+outofAFRpie & theme(legend.position = "right")
combined + plot_layout(ncol=2, nrow=1, guides = "collect", axis_titles = "collect", tag_level="keep")

TE_africaoutofafrica_anno_counted <- rbind(TE_outofafrica_anno_counted,TE_africa_anno_counted)

ggplot(TE_africaoutofafrica_anno_counted, mapping=aes(x=pop, y=NUMBER,
                                        fill=fct_relevel(ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  theme_classic()

# PCA ---------------------------------------------------------------------

## define variables and plink output files

pca <- read_table2(file.path(sys_dir,"PanTE_human/pca/TE_pca.eigenvec"), col_names = FALSE)
eigenval <- scan(file.path(sys_dir,"PanTE_human/pca/TE_pca.eigenval"))

AFR <- c("HG01891","HG02257","HG02486","HG02559","HG02572","HG03516","HG02622","HG02630","HG02717","HG02886","HG03453","HG03540","HG03579","HG02109","HG02145","HG02723","HG02818","HG03486","NA18906","NA20129","NA21309","HG02055","HG03098")
outofAFR <- c("HG01123","HG01258","HG01358","HG01361","HG00735","HG00741","HG01071","HG01106","HG01175","HG01928","HG01952","HG01978","HG02148","HG00733","HG01109","HG01243","HG00438","HG00621","HG00673","HG02080","HG03492")
pop_info <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/hprc_year1_sample_metadata.txt"), header = F)
colnames(pop_info) <- c("ind","subpop","pop")
subpop_vectors <- list()
unique_subpops <- unique(pop_info$subpop)
for (subpop in unique_subpops) {
  subpop_inds <- pop_info[pop_info$subpop == subpop, "ind"]
  subpop_vectors[[subpop]] <- subpop_inds}

## housekeeping for pca variables

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# spp <- rep(NA, length(pca$ind))
# spp[grep("", pca$ind)] <- "human"
pop <- rep(NA, length(pca$ind))
pop[grep("CHM13", pca$ind)] <- "CHM13"
for (value in AFR) {pop[grep(value, pca$ind)] <- "Africa"}
for (value in outofAFR) {pop[grep(value, pca$ind)] <- "outofAfrica"}
# spp_pop <- paste0(spp, "_", pop)
subpop <- rep(NA, length(pca$ind))
subpop[grep("CHM13", pca$ind)] <- "CHM13"
for (value in  subpop_vectors[["CLM"]] ) {subpop[grep(value, pca$ind)] <- "CLM"}
for (value in  subpop_vectors[["ACB"]] ) {subpop[grep(value, pca$ind)] <- "ACB"}
for (value in  subpop_vectors[["GWD"]] ) {subpop[grep(value, pca$ind)] <- "GWD"}
for (value in  subpop_vectors[["ESN"]] ) {subpop[grep(value, pca$ind)] <- "ESN"}
for (value in  subpop_vectors[["CHS"]] ) {subpop[grep(value, pca$ind)] <- "CHS"}
for (value in  subpop_vectors[["PUR"]] ) {subpop[grep(value, pca$ind)] <- "PUR"}
for (value in  subpop_vectors[["PEL"]] ) {subpop[grep(value, pca$ind)] <- "PEL"}
for (value in  subpop_vectors[["MSL"]] ) {subpop[grep(value, pca$ind)] <- "MSL"}
for (value in  subpop_vectors[["KHV"]] ) {subpop[grep(value, pca$ind)] <- "KHV"}
for (value in  subpop_vectors[["PJL"]] ) {subpop[grep(value, pca$ind)] <- "PJL"}
for (value in  subpop_vectors[["YRI"]] ) {subpop[grep(value, pca$ind)] <- "YRI"}
for (value in  subpop_vectors[["ASW"]] ) {subpop[grep(value, pca$ind)] <- "ASW"}
for (value in  subpop_vectors[["KENYA"]] ) {subpop[grep(value, pca$ind)] <- "KENYA"}

pca <- as.tibble(data.frame(pca, pop, subpop)) # removed spp and spp_pop because all ind are humans

## plot the % of explained variance and the PCA

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)

pop_order <- fct_relevel(pca$pop,"Africa","outofAfrica","CHM13")

ggplot(pca, aes(PC1, PC2, col = pop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("#00539CFF","#EEA47FFF","#A2A2A1FF"))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  labs(color="Population",shape="Population")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pca_nochm13 <- pca %>% filter(pca$ind != "CHM13")
pop_order <- fct_relevel(pca_nochm13$pop,"Africa","outofAfrica")
ggplot(pca_nochm13, aes(PC1, PC2, col = pop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("#00539CFF","#EEA47FFF"))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  labs(color="Population",shape="Population")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_pca.png"), plot=ggplot2::last_plot())

subpop_order <- fct_relevel(pca$subpop,"ACB","GWD","MSL","ESN","YRI","ASW","KENYA",
                            "CHS","PUR","CLM","PEL","KHV","PJL",
                            "CHM13")
pop_order <- fct_relevel(pca$pop,"Africa","outofAfrica","CHM13")

ggplot(pca, aes(PC1, PC2, col = subpop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_color_viridis(discrete = TRUE, option = "C")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  labs(color="Subpopulation",shape="Population")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_pca.png"), plot=ggplot2::last_plot())

# Frequency based selection scans -----------------------------------------

Fst <- read.csv(file.path(sys_dir,"PanTE_human/selection_scans/fst/out.weir.fst"), header = T,sep = "")
Fst <- Fst %>% filter(Fst$WEIR_AND_COCKERHAM_FST >= 0)
Fst$CHROM <- gsub("chr","",Fst$CHROM)
Fst$CHROM <- as.numeric(Fst$CHROM)
Fstsubset <- Fst[complete.cases(Fst),]
TE <- c(1:(nrow(Fstsubset)))
mydf <- data.frame(TE,Fstsubset)
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")
mydf_99 <- quantile(mydf$WEIR_AND_COCKERHAM_FST, 0.99)
mydf_sig <- mydf[mydf$WEIR_AND_COCKERHAM_FST > mydf_99,] 
manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="TE",logp=FALSE,ylab="Fst", cex=2, xaxt="n",
          col=c(gray_chr_colors))
abline(h = mydf_99, col = "#EE353E", lty = 2)

# Haplotype based selection scans -----------------------------------------

## rehh - define data files

rehh_wrapper <- function(vcf_file, chr_name) {
  result <- data2haplohh(hap_file = vcf_file,
                         chr.name = chr_name,
                         polarize_vcf = F,
                         verbose = FALSE,
                         recode.allele = FALSE,
                         vcf_reader = "vcfR")
  scan <- scan_hh(result)
  return(scan)
  # assign(paste0("AFR_", chr_name), result, envir = .GlobalEnv)
  cat("Converting and scanning data:", vcf_file, "\n")
}

## outofAFR

directory_path <- file.path(sys_dir,"PanTE_human/selection_scans/ihs_xpehh/outofAFR")
file_pattern <- "outofAFR_1alt_phased_"
file_list <- list.files(directory_path, pattern = file_pattern, full.names = TRUE)
chr_id = as.factor(c("chr1", "chr10" ,"chr11" ,"chr12" ,"chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", "chr18", "chr19" ,"chr2", "chr20","chr21" ,"chr22" ,"chr3" ,"chr4", "chr5" ,"chr6" ,"chr7" ,"chr8" ,"chr9"))
rehh_outofAFR <- data.frame() 

for (vcf_file in file_list) {
  for (chr_level in levels(chr_id)) {
    chr_name <- chr_level
    chr_name <- as.character(chr_level)
    if (grepl(paste0("phased_", chr_name, "\\.vcf"), vcf_file)) {
      # results_list_outofAFR[[paste0("outofAFR_", chr_name)]] <- get(paste0("outofAFR_", chr_name), envir = .GlobalEnv)
      scan <- rehh_wrapper(vcf_file, chr_name)
      rehh_outofAFR <- rbind(rehh_outofAFR, scan)  }}}

## AFR

directory_path <- file.path(sys_dir,"PanTE_human/selection_scans/ihs_xpehh/AFR")
file_pattern <- "AFR_1alt_phased_"
file_list <- list.files(directory_path, pattern = file_pattern, full.names = TRUE)
chr_id = as.factor(c("chr1", "chr10" ,"chr11" ,"chr12" ,"chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", "chr18", "chr19" ,"chr2", "chr20","chr21" ,"chr22" ,"chr3" ,"chr4", "chr5" ,"chr6" ,"chr7" ,"chr8" ,"chr9"))
rehh_AFR <- data.frame() 

for (vcf_file in file_list) {
  for (chr_level in levels(chr_id)) {
    chr_name <- chr_level
    chr_name <- as.character(chr_level)
    if (grepl(paste0("phased_", chr_name, "\\.vcf"), vcf_file)) {
      scan <- rehh_wrapper(vcf_file, chr_name) 
      # results_list_AFR[[paste0("AFR_", chr_name)]] <- get(paste0("AFR_", chr_name), envir = .GlobalEnv)
      rehh_AFR <- rbind(rehh_AFR, scan)  }}}

rehh_AFR
rehh_outofAFR

rehh_AFR$CHR <- gsub("chr1([^0-9]|$)", "chr01\\1", rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr2([^0-9]|$)", "chr02\\1", rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr3","chr03",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr4","chr04",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr5","chr05",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr6","chr06",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr7","chr07",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr8","chr08",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr9","chr09",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr","",rehh_AFR$CHR)
rehh_AFR$CHR <- as.numeric(rehh_AFR$CHR)

rehh_outofAFR$CHR <- gsub("chr1([^0-9]|$)", "chr01\\1", rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr2([^0-9]|$)", "chr02\\1", rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr3","chr03",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr4","chr04",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr5","chr05",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr6","chr06",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr7","chr07",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr8","chr08",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr9","chr09",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr","",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- as.numeric(rehh_outofAFR$CHR)

chr_colors <- brewer.pal(n = 8, name = "Dark2")

## ihs - calculate ihs and manhattan plots

ihs_AFR <- ihh2ihs(rehh_AFR, min_maf = -1)
palette(chr_colors)
ihs_AFR_2 <- ihs_AFR$ihs %>% drop_na
ihs_AFR_2$IHS_abs <- abs(ihs_AFR_2$IHS)
per_99 <- quantile(ihs_AFR_2$IHS, 0.99)
per_01 <- quantile(ihs_AFR_2$IHS, 0.01)
manhattanplot(ihs_AFR,
              cex = 1.5,
              threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20","21","22"))
freqbinplot(ihs_AFR)
# manhattanplot(ihs_AFR,
#               pval = TRUE,
#               threshold = 2)

ihs_outofAFR <- ihh2ihs(rehh_outofAFR, min_maf = -1) 
palette(chr_colors)
ihs_outofAFR_2 <- ihs_outofAFR$ihs %>% drop_na
ihs_outofAFR_2$IHS_abs <- abs(ihs_outofAFR_2$IHS)
per_99 <- quantile(ihs_outofAFR_2$IHS, 0.99)
per_01 <- quantile(ihs_outofAFR_2$IHS, 0.01)
manhattanplot(ihs_outofAFR,
              cex = 1.5,
              threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20","21","22"))
freqbinplot(ihs_outofAFR)
# manhattanplot(ihs_outofAFR,
#               pval = TRUE,
#               threshold = 2)

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

ihs_AFR_outofAFR <- rbind(ihs_outofAFR_2,ihs_AFR_2)
ihs_AFR_outofAFR$TE <- seq(1,nrow(ihs_AFR_outofAFR))
# ihs_AFR_outofAFR$IHS_abs <- abs(ihs_AFR_outofAFR$IHS) 

manhattan(ihs_AFR_outofAFR,chr="CHR",bp="POSITION",p="IHS_abs",snp="TE",logp=FALSE,ylab="iHS", cex=2, xaxt="n",
          col=c(gray_chr_colors),suggestiveline=F,genomewideline=F)
per_99 <- quantile(ihs_AFR_2$IHS_abs, 0.99)
abline(h = per_99, col = "#EE353E", lty = 2)
per_99 <- quantile(ihs_outofAFR_2$IHS_abs, 0.99)
abline(h = per_99, col = "steelblue", lty = 2)

## ihs - extract candidate TE loci with highest 99% absolute ihs

ihs_outofAFR <- ihs_outofAFR$ihs %>% drop_na
ihs_outofAFR$IHS_abs <- abs(ihs_outofAFR$IHS)
per_99 <- quantile(ihs_outofAFR$IHS_abs, 0.99)
ihs_outofAFR_per_99 <- ihs_outofAFR[ihs_outofAFR$IHS_abs > per_99,]
ihs_outofAFR_candidate_TE <- ihs_outofAFR_per_99
# per_99 <- quantile(ihs_outofAFR$IHS, 0.99)
# per_01 <- quantile(ihs_outofAFR$IHS, 0.01)
# ihs_outofAFR_per_01 <- ihs_outofAFR[ihs_outofAFR$IHS < per_01,]
# ihs_outofAFR_per_99 <- ihs_outofAFR[ihs_outofAFR$IHS > per_99,]
# ihs_outofAFR_candidate_TE <- ihs_outofAFR %>% 
#   filter(LOGPVALUE > 2 | LOGPVALUE < -2) %>% 
#   dplyr::select(c("CHR","POSITION"))

ihs_AFR <- ihs_AFR$ihs %>% drop_na
ihs_AFR$IHS_abs <- abs(ihs_AFR$IHS)
per_99 <- quantile(ihs_AFR$IHS_abs, 0.99)
ihs_AFR_per_99 <- ihs_AFR[ihs_AFR$IHS_abs > per_99,]
ihs_AFR_candidate_TE <- ihs_AFR_per_99
# per_99 <- quantile(ihs_AFR$IHS, 0.99)
# per_01 <- quantile(ihs_AFR$IHS, 0.01)
# ihs_AFR_per_01 <- ihs_AFR[ihs_AFR$IHS < per_01,]
# ihs_AFR_per_99 <- ihs_AFR[ihs_AFR$IHS > per_99,]
# ihs_AFR_candidate_TE <- ihs_AFR %>% 
#   filter(LOGPVALUE > 2 | LOGPVALUE < -2) %>% 
#   dplyr::select(c("CHR","POSITION"))

ihs_AFR_outofAFR_candidate_TE <- rbind(ihs_AFR_candidate_TE,ihs_outofAFR_candidate_TE)
ihs_AFR_outofAFR_candidate_TE <- ihs_AFR_outofAFR_candidate_TE[order(ihs_AFR_outofAFR_candidate_TE[,1],ihs_AFR_outofAFR_candidate_TE[,2]),]
ihs_AFR_outofAFR_candidate_TE$POSITION <- as.numeric(ihs_AFR_outofAFR_candidate_TE$POSITION)
write.table(ihs_AFR_outofAFR_candidate_TE, file = file.path(sys_dir,'PanTE_human/ihs_AFR_outofAFR_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

## xpehh - calculate xpehh and manhattan plots

xpehh <- ies2xpehh(scan_pop1 =  rehh_AFR,
                   scan_pop2 =  rehh_outofAFR,
                   popname1 = "AFR",
                   popname2 = "outofAFR") %>% drop_na()

manhattanplot(xpehh,
              cex = 1.5,
              # threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20"))

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

xpehh$TE <- seq(1,nrow(xpehh))

manhattan(xpehh,chr="CHR",bp="POSITION",p="XPEHH_AFR_outofAFR",snp="TE",logp=FALSE,ylab="xp-EHH", cex=2, xaxt="n",
          col=c(gray_chr_colors))
abline(h = per_99, col = "#EE353E", lty = 2)

## xpehh - extract candidate TE loci with highest 99% absolute ihs

per_99 <- quantile(xpehh$XPEHH_AFR_outofAFR, 0.99)
per_01 <- quantile(xpehh$XPEHH_AFR_outofAFR, 0.01)
xpehh_per_01 <- xpehh[xpehh$XPEHH_AFR_outofAFR < per_01,]
xpehh_per_99 <- xpehh[xpehh$XPEHH_AFR_outofAFR > per_99,]
xpehh_per_0199 <- rbind(xpehh_per_99,xpehh_per_01)
xpehh_per_0199_new_candidate_TE <- xpehh_per_0199 %>% dplyr::select(c("CHR","POSITION"))
xpehh_per_0199_new_candidate_TE <- xpehh_per_0199_new_candidate_TE[order(xpehh_per_0199_new_candidate_TE[,1],xpehh_per_0199_new_candidate_TE[,2]),]
xpehh_per_0199_new_candidate_TE$POSITION <- as.numeric(xpehh_per_0199_new_candidate_TE$POSITION)
write.table(xpehh_per_0199_new_candidate_TE, file = file.path(sys_dir,'PanTE_human/xpehh_per_0199_new_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

# Ohana demography based selection scans ----------------------------------

ohana <- read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/scansel.txt"), header=T)
pos <-  read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/TE.map"), header=F) %>% dplyr::select(c(4))
all_pos <-  read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/input_pos_ohana.txt"), header=F)
colnames(all_pos) <- c("chr","pos")
ohana_pos <- cbind(ohana,pos) %>% dplyr::select(c(-1)) %>% data.frame()
colnames(ohana_pos) <- c("global.lle","local.lle","lle.ratio","pos")
ohana_postions <- merge(ohana_pos,all_pos)

percentile_99 <- quantile(ohana_postions$lle.ratio, 0.99)
ohana_pos_sig <- ohana_postions %>% filter(ohana_postions$lle.ratio > percentile_99)
ohana_pos_sig_pos <- ohana_pos_sig %>% dplyr::select("chr","pos")
write.table(ohana_pos_sig_pos, file = file.path(sys_dir,'PanTE_human/ohana_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

sig <-  ohana_postions %>% filter(ohana_postions$lle.ratio > percentile_99)
nrow(sig)

ohana_postions$chr <- gsub("chr1","1",ohana_postions$chr)
ohana_postions$chr <- gsub("chr2","2",ohana_postions$chr)
ohana_postions$chr <- gsub("chr(1[1-9])([^0-9]|$)", "\\1\\2", ohana_postions$chr)
ohana_postions$chr <- gsub("chr3","3",ohana_postions$chr)
ohana_postions$chr <- gsub("chr4","4",ohana_postions$chr)
ohana_postions$chr <- gsub("chr5","5",ohana_postions$chr)
ohana_postions$chr <- gsub("chr6","6",ohana_postions$chr)
ohana_postions$chr <- gsub("chr7","7",ohana_postions$chr)
ohana_postions$chr <- gsub("chr8","8",ohana_postions$chr)
ohana_postions$chr <- gsub("chr9","9",ohana_postions$chr)
ohana_postions$chr <- gsub("chr20","20",ohana_postions$chr)
ohana_postions$chr <- gsub("chr21","21",ohana_postions$chr)
ohana_postions$chr <- gsub("chr22","22",ohana_postions$chr)
ohana_postions$chr <- gsub("chrX","22",ohana_postions$chr)
ohana_postions$chr <- gsub("chrY","23",ohana_postions$chr)
unique(ohana_postions$chr)
ohana_postions$chr <- as.numeric(ohana_postions$chr)
ohana_postions$TE <- seq(1,nrow(ohana_postions))
ohana_postions$lle.ratio_log <- log10(ohana_postions$lle.ratio)
percentile_99 <- quantile(ohana_postions$lle.ratio_log, 0.99)

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

manhattan(ohana_postions,chr="chr",bp="pos",p="lle.ratio_log",snp="TE",logp=F,ylab="ohana (log scale)", cex=2, xaxt="n",
          col=c(gray_chr_colors),suggestiveline=F,genomewideline=F)
abline(h = (percentile_99), col = "#EE353E", lty = 2)

# Candidate loci - All selection scans ------------------------------------

## extract and merge candidate TE loci from all selection scans 

candidate_fst <- candidate_TE_div_Fst %>% dplyr::select("CHR","POS")
candidate_ihs <- ihs_AFR_outofAFR_candidate_TE %>% dplyr::select("CHR","POSITION")
candidate_xpehh <- xpehh_per_0199_new_candidate_TE %>% dplyr::select("CHR","POSITION")
candidate_ohana <- ohana_pos_sig_pos %>% dplyr::select("chr","pos")
colnames(candidate_fst) <- c("chr","pos")
colnames(candidate_ihs) <- c("chr","pos")
colnames(candidate_xpehh) <- c("chr","pos")
colnames(candidate_ohana) <- c("chr","pos")
candidate_fst$selscan <- "fst"
candidate_ihs$selscan <- "ihs"
candidate_xpehh$selscan <- "xpehh"
candidate_ohana$selscan <- "ohana"

candi_all <- rbind(candidate_fst,candidate_ihs,candidate_xpehh,candidate_ohana) %>% data.frame() 
candi_all$chr <- gsub("chr1","1",candi_all$chr)
candi_all$chr <- gsub("chr2","2",candi_all$chr)
candi_all$chr <- gsub("chr(1[1-9])([^0-9]|$)", "\\1\\2", candi_all$chr)
candi_all$chr <- gsub("chr3","3",candi_all$chr)
candi_all$chr <- gsub("chr4","4",candi_all$chr)
candi_all$chr <- gsub("chr5","5",candi_all$chr)
candi_all$chr <- gsub("chr6","6",candi_all$chr)
candi_all$chr <- gsub("chr7","7",candi_all$chr)
candi_all$chr <- gsub("chr8","8",candi_all$chr)
candi_all$chr <- gsub("chr9","9",candi_all$chr)
candi_all$chr <- gsub("chr20","20",candi_all$chr)
candi_all$chr <- gsub("chr21","21",candi_all$chr)
candi_all$chr <- gsub("chr22","22",candi_all$chr)
candi_all <- candi_all %>% mutate(chr = paste0("chr", chr))
candi_all <- candi_all[order(candi_all[,1],candi_all[,2]),]
rownames(candi_all) <- 1:nrow(candi_all)

duplicated_entries <- duplicated(candi_all[c("chr", "pos")]) | duplicated(candi_all[c("chr", "pos")], fromLast = TRUE)
candi_all_2 <- candi_all[!duplicated_entries, ]

write_xlsx(candi_all_2, file.path(sys_dir,'PanTE_human/candi_all.xlsx'))
write.table(candi_all_2, file.path(sys_dir,'PanTE_human/candi_all.bed'),sep = '\t', row.names = F,quote=F)

## merge with TE_div_Fst_filtered and TE_FULL_INFO to get full info on candidate TE loci

colnames(TE_div_Fst_filtered)
colnames(candi_all_2) <- c("CHR","POS","selscan")

candi_all_2$POS <- as.numeric(candi_all_2$POS)

bin          <- 10000

df2          <- candi_all_2 %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_selscan = list(selscan))
df2[,3:4] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_selscan","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_selscan"))
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           sum_selscan = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
candi_all_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

colnames(TE_FULL_INFO)
colnames(candi_all_binned)
colnames(TE_size_distTSS_div_AO)
colnames(candi_all)

candi_div_Fst <- merge(TE_div_Fst_filtered,candi_all)
candi_size_distTSS_div_AO <- merge(TE_size_distTSS_div_AO,candi_all)
candi_FULL_INFO <- merge(TE_FULL_INFO,candi_all_binned)

summary(candi_div_Fst$FREQ_ALT)
summary(candi_div_Fst$SIZE)
table(candi_div_Fst$anno_order)
table(candi_div_Fst$selscan)

summary(candi_size_distTSS_div_AO$FREQ_ALT)
summary(candi_size_distTSS_div_AO$SIZE)
table(candi_size_distTSS_div_AO$ANNO)
table(candi_size_distTSS_div_AO$selscan)

candi_FULL_INFO_HR <- candi_FULL_INFO %>% filter(mean_recrate > 5)

summary(candi_FULL_INFO_HR$mean_recrate)
summary(candi_FULL_INFO_HR$mean_disttss)








# Supp Tables -------------------------------------------------------------

Rdataframes_CT <- data.frame(
    Sheet = c("Sheet1", "Sheet2", "Sheet3", "Sheet4", "Sheet5", "Sheet6", 
              "Sheet7", "Sheet8", "Sheet9", "Sheet10", "Sheet11", "Sheet12", 
              "Sheet13", "Sheet14", "Sheet15", "Sheet16", "Sheet17", "Sheet18"),
    DataFrameName = c("global_freq", "global_freq_annotation", "genomic_location", 
                      "TE_dataframe", "TE_divergence_fst", "TE_FULL_INFO", 
                      "TE_size_chromosome_length_counted", "TE_size_chromosome_length", "TE_counts_per_chromosome", 
                      "background_recrate", "major_singletons_shared_per_individual", "TE_counts_per_individual", 
                      "ihs", "xpehh", "ohana", "candidate_TE_ihs", "candidate_TE_xpehh", "candidate_TE_all"))

write_xlsx(list(
                ContentTable = Rdataframes_CT,
                Sheet1 = super_freq_7,
                Sheet2 = super_freq_anno_counted,
                Sheet3 = TE_loc_c,
                Sheet4 = TE_dataframe,
                Sheet5 = TE_div_Fst_filtered,
                Sheet6 = TE_FULL_INFO,
                Sheet7 = size_length_counted,
                Sheet8 = size_length,
                Sheet9 = biTE_count_2,
                Sheet10 = background_rec,
                Sheet11 = maj_sin_shar_tidy,
                Sheet12 = TE_indv,
                Sheet13 = ihs_AFR_outofAFR,
                Sheet14 = xpehh,
                Sheet15 = ohana_postions,
                Sheet16 = ihs_AFR_outofAFR_candidate_TE,
                Sheet17 = xpehh_per_0199,
                Sheet18 = candi_all),
           file.path(sys_dir,"PanTE_human/scripts/supp_data_Rdataframes.xlsx"))







