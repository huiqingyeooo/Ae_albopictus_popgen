# Identity by descent to identify close kins

library(SNPRelate)
library(xlsx)

# Read genotypes
bed.fn <- "albo_prune.bed"
fam.fn <- "albo_prune.fam"
bim.fn <- "albo_prune.bim"

# Convert to gds format
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "albo_ibd.gds")
albo <- snpgdsOpen("albo_ibd.gds")
snpgdsSummary("albo_ibd.gds")

# Identity by descent analysis
albo_ibd <- snpgdsIBDMLE(albo, autosome.only = FALSE, kinship = TRUE)
#albo_ibd <- snpgdsIBDMLE(albo, autosome.only = FALSE, kinship = TRUE, num.thread=2)

# Output coefficients
albo_coefficients <- snpgdsIBDSelection(albo_ibd)

# Plot
tiff("coefficients.tiff", width = 15, height = 10, units = 'in', res=100)
plot(albo_coefficients$k0, albo_coefficients$k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", main="albo_coefficients")
lines(c(0,1), c(1,0), col="red", lty=2)
dev.off()

#k0 - refers to the probability of two diploid individuals sharing 0 alleles that are IBD
#k1 - prob. of two diploid individuals sharing 1 allele that are IBD

write_xlsx(albo_coefficients, "./albo_coefficients.xlsx")

