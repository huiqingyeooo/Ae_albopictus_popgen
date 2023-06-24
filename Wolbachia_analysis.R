# Wolbachia infection analysis

library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)

#####################################################################
# Test if Wolbachia loads of Aedes albopictus differ across habitats
#####################################################################

qPCR <- read.csv("compiled_qPCR.csv", header=T)
qPCR$habitat<-as.factor(qPCR$habitat)

# Remove points with >0.5 difference
qPCR <- subset(qPCR, !Remarks=="difference_0.5")
qPCR <- subset(qPCR, !Remarks=="exclude")

# Identify and remove outliers
qPCR$wolB_density #outlier datapoint at row 37
wolB <- qPCR[-c(37),]

# Statistical tests
kruskal.test(wolB_density~habitat, data=wolB) #compare between > 2 groups
pairwise.wilcox.test(wolB$wolB_density, wolB$habitat, p.adjust.method = "bonferroni") #pairwise comparison

# Quick summary of samples
wolB%>%
  group_by(habitat) %>%
  tally()

# Boxplot
Bplot <- wolB%>%
  ggplot(aes(x=habitat, y=wolB_density)) +
  geom_jitter(aes(color=habitat),size=3,alpha=0.5,width=0.15)+
  geom_boxplot(alpha=0) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_colour_manual(values=c("#99D98C","#FFC100","#168AAD"))+
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  xlab("Habitat type") + ylab("wAlbB density")+
  ylim(0,33)

# Perform pairwise comparisons and add signficance values to plot
compare_means(wolB_density ~ habitat,  data = wolB)

# Visualize: Specify the comparisons you want
my_comparisons <- list(c("F", "P"), c("P", "U"))
Bplot+ 
  stat_compare_means(comparisons = my_comparisons, label.y=c(30,35))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

# Plot for manuscript
tiff(filename="wolB_boxplot_color.tiff",width=3,height=3.5,res=300,units="in")
Bplot+ 
  stat_compare_means(comparisons = my_comparisons, label.y=c(27,29))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 33, label.x=0.5)     # Add global p-value
dev.off()

#################################################
# Test if Wolbachia A and B loads are correlated
#################################################
# Standardised Major Axis Estimation
library(smatr)

qPCR <- read.csv("compiled_qPCR.csv", header=T)

# Remove points with >0.5 difference
qPCR <- subset(qPCR, !Remarks=="difference_0.5")
qPCR <- subset(qPCR, !Remarks=="exclude")

# Identify and remove outliers
qPCR$wolB_density #outlier datapoint at row 37
qPCR$wolA_density #outlier datapoint at row 16
wolAB <- qPCR[-c(16,37),]

# Transform datapoints with log
wolAB$logWolA<-log10(wolAB$wolA_density+0.01)
wolAB$logWolB<-log10(wolAB$wolB_density+0.01)

# SMA model
model <- sma(logWolB~logWolA, data=wolAB,test="elevation")
summary(model)

# SMA plot, with points coloured by habitat type
tiff(filename="wolA_vs_wolB_sma.tiff",width=6,height=5,res=300,units="in")
plot(model,ylab="Log(Wolbachia B density)",xlab="Log(Wolbachia A density)",
     col="#616161",lwd=2,pch=20)
F<-wolAB %>% filter(habitat=="F")
P<-wolAB %>% filter(habitat=="P")
U<-wolAB %>% filter(habitat=="U")
points(F$logWolB~F$logWolA, col=c("#99D98C"),pch=16,cex=1.2)
points(P$logWolB~P$logWolA, col=c("#FFC100"),pch=16,cex=1.2)
points(U$logWolB~U$logWolA, col=c("#168AAD"),pch=16,cex=1.2)
abline(model, col="#616161",lwd=2)
dev.off()