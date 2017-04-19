
setwd("/xchip/cga_home/margolis/retainedIntron/Snyder_ipi/kma")
#getwd()

#Sys.setenv(R_LIBS_USER = "/xchip/cga_home/asmart/KMA")
#Sys.getenv("R_LIBS_USER")
#.libPaths("/xchip/cga_home/asmart/KMA")
#.libPaths()

#install.packages("git2r", repos="http://cran.us.r-project.org")
#required_packages <- c("devtools", "data.table", "reshape2", "dplyr")
#install.packages(required_packages, repos='http://cran.us.r-project.org')

#devtools::install_github("pachterlab/kma")

library("kma", lib.loc="/xchip/cga_home/margolis/Packages/KMA")
library(methods)
print("library loaded")

#system.file("pre-process", package="kma")

print("loading express files")
xprsDir = file.path("/xchip/cga_home/margolis/retainedIntron/Snyder_ipi")
xprs_fnames = list.files(xprsDir, pattern="results_3.xprs", full = TRUE, recursive = TRUE)
print("express filenames")
xprs_fnames

print("loading decoder")
decoder <- read.table("/xchip/cga_home/margolis/retainedIntron/Snyder_ipi/kma/decoder.txt", header=T, sep='\t', stringsAsFactors=FALSE)
decoder$group.ID <- factor(decoder$group.ID, levels = c("LB", "NB"))
print("set factor levels: LB, NB")
decoder

sample_names = decoder$sample.ID
condition_names = decoder$group.ID

print("reading express files")
xprs <- read_express(xprs_fnames, sample_names, condition_names)
names(xprs)
save.image()

#load(".RData")
#print("loaded RData")

intron_to_trans <- data.table::fread(file.path("/xchip/cga_home/asmart/KMA/160705_output", "intron_to_transcripts_v2.txt"), data.table = FALSE)
head(intron_to_trans)

ir <- newIntronRetention(xprs$tpm, intron_to_trans, xprs$condition, xprs$uniq_counts)
print(ir)
head(ir$flat)

ir <- ir %>%
    filter_low_tpm(1) %>%
    filter_perfect_psi() %>%
    filter_low_frags(3)

colnames(ir$flat)
save.image()

##need to generate zero coverage file using python
print("loading zero coverage files")
zc_dir = file.path("/xchip/cga_home/margolis/retainedIntron/Snyder_ipi")
zc_fnames = list.files(zc_dir, pattern="zero_coverage_3.txt", full = TRUE, recursive = TRUE)
print("zc_fnames:")
zc_fnames

decoder <- read.table("/xchip/cga_home/margolis/retainedIntron/Snyder_ipi/kma/decoder.txt", header=T, sep='\t', stringsAsFactors=FALSE)

zc_samples = decoder$sample.ID
zc_conditions = decoder$group.ID
print("loaded decoder")

all_zc <- get_batch_intron_zc(zc_fnames, zc_samples, zc_conditions)
print("all zc:")
head(all_zc)
save.image()

ir <- summarize_zero_coverage(ir, all_zc)
print("summarized zero coverage")

#set.seed(42)
ir_test <- retention_test(ir)
## aggregating filters
## joining filtered data set
## estimating the null for each condition
## computing p-values
print("ir test:")
head(ir_test)


ir_test %>%
    filter(qvalue <= 0.10) %>%
    select(-c(pvalue))


save.image()
write.csv(ir$flat, file="160922_Snyder_ipi_KMA_IR_flat_filtered.csv")
write.csv(ir_test, file="160922_Snyder_ipi__KMA_IR_test.csv")
