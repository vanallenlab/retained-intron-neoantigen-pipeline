library(sleuth)

library(devtools)
dev_mode()
print("Dev_mode() on")

setwd("XXXXXX")  # Set your working directory here

library("kma", lib.loc="XXXXX")  # lib.loc = wherever you downloaded KMA
library(methods)
print("library loaded")

kallistoDir = file.path("XXXXX")  # Directory where your kallisto output files live 
kallisto_fnames = list.files(kallistoDir, pattern="abundance.tsv", full = TRUE, recursive = TRUE)
print("kallisto filenames")
kallisto_fnames

print("loading decoder")
decoder <- read.table("decoder_all.txt", header=T, sep='\t', stringsAsFactors=FALSE)
decoder$group.ID <- factor(decoder$group.ID, levels = c("XXXXX","XXXXX","XXXXX")) # Insert your appropriate group levels here
print("set factor levels")
decoder

sample_names = decoder$sample.ID
condition_names = decoder$group.ID

print("reading kallisto files")
kallisto <- read_kallisto_kma(kallisto_fnames, sample_names, condition_names, cols = c('tpm', 'est_counts'), substring_replace = '\\|\\S*$')  # NOTE: Decoder must be in same order as list.files for this to work properly!!!!!! 
names(kallisto)
save.image()

intron_to_trans <- data.table::fread(file.path("XXXXXX", "intron_to_transcripts.txt"), data.table = FALSE)  # Insert path to the intron_to_transcripts.txt file that you created in the KMA pre-processing step
head(intron_to_trans)

ir <- newIntronRetention(kallisto$tpm, intron_to_trans, kallisto$condition, kallisto$est_counts)
print(ir)
head(ir$flat)

n_filter <- 0.25 # Set chosen sample percentage here
ir <- filter_low_tpm(ir, 1, n_filter) # Only include transcripts that are expressed at a tpm >= 1 in > n_filter % of samples
ir <- filter_perfect_psi(ir) # Remove transcripts in which every sample has 0 or 1 retention after rounding (artifact)
ir <- filter_low_frags(ir, 5, n_filter) # Only include transcripts that have at least 5 unique counts in > n_filter % of samples
colnames(ir$flat)

write.csv(ir$flat, file="Aggregate_KMA-kallisto_flat_filtered.csv")  # This file is output for the retained intron neoantigen downstream pipeline


print("aggregating filters")
agg_ir <- aggregate_filters(ir)
agg_ir <- dplyr::filter(agg_ir, f_all)

valid_intron <- unique(agg_ir$intron)
filtered_retention <- ir$retention[valid_intron, ]


set.seed(42)
ir_test <- retention_test(ir)
print("ir test:")
head(ir_test)

ir_test %>%
    filter(qvalue <= 0.10) %>%
    select(-c(pvalue))

write.csv(ir_test, file="Aggregate_KMA-kallisto_test.csv")
