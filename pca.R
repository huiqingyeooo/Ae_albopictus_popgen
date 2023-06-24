# Principal Component Analysis
rm(list = ls())

library(SNPRelate)
library(writexl)
library(readxl)
library(ggplot2)
library(ggrepel)

albo <- "./albo_plink.vcf"
snpgdsVCF2GDS(albo, "albo.gds", method = "copy.num.of.ref")
snpgdsSummary("albo.gds")

albo.gdsfile <- snpgdsOpen("albo.gds")
albo.pca <- snpgdsPCA(albo.gdsfile, autosome.only = FALSE)

# % of variation accounted for by top principal components
pc.percent <-albo.pca$varprop*100
head(round(pc.percent, 2)) #Take note of values for PC1 and PC2

# Calculate correlations between eigenvectors and SNP genotypes
corr <- snpgdsPCACorr(albo.pca,albo.gdsfile,with.id=TRUE)
hist(corr$snpcorr[1,])

scaffold <- read.gdsn(index.gdsn(albo.gdsfile, "snp.chromosome"))
length(scaffold) #27037
plot(abs(corr$snpcorr[2,]), xlab="SNP Index", ylab="PC 2")#,col=scaffold)

# Make a data frame to be extracted
tab <- data.frame (sample.id = albo.pca$sample.id,
                   EV1 = albo.pca$eigenvect[,1],EV2 = albo.pca$eigenvect[,2],EV3 = albo.pca$eigenvect[,3],
                   EV4 = albo.pca$eigenvect[,4],EV5 = albo.pca$eigenvect[,5],EV6 = albo.pca$eigenvect[,6],
                   EV7 = albo.pca$eigenvect[,7],EV8 = albo.pca$eigenvect[,8],EV9 = albo.pca$eigenvect[,9],
                   EV10 = albo.pca$eigenvect[,10],stringsAsFactors = FALSE)

# Extract factors for external plotting
write_xlsx(tab, "albo_pca.xlsx")

# After saving excel file, insert column with the heading population

# Plot PCA
pca <- read_xlsx("albo_pca.xlsx", 1)
ggdata <- data.frame(pca)
str(ggdata)

habitat_colour<- c("#99D98C","#FFC100","#168AAD")

pca.plot <- ggplot(ggdata) +
  geom_point(aes(x=EV1, y=EV2, color=habitat), size=5,show.legend = FALSE) +
  scale_color_manual(values = habitat_colour, name = "habitat") +
  theme_bw(base_size=12) +
  labs(x ="PC1 (1.38%)", y = "PC2 (1.27%)") +  # remember to change this section!
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        axis.text.x=element_text (size=20),
        axis.text.y=element_text (size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20, face="bold"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        plot.title = element_text(size=25))

tiff("albo_habitat.tiff", width = 15, height = 10, units = 'in', res=300)
pca.plot
dev.off()
