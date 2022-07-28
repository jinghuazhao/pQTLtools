# MR analysis

# load libraries and helper functions
source("workflow/scripts/MR_functions.R")

# read & reformat data
data_MR <- fread(snakemake@input[["data_MR"]])

setnames(data_MR, c("rsID", "beta_Prot", "se_Prot", "beta_HF", "se_HF"),
         c("SNP", "beta.x", "se.x", "beta.y", "se.y"))

snplist <- read_lines(snakemake@input[["snplist"]])

ldrho <- fread(snakemake@input[["ldrho"]],
               col.names = snplist)  %>%
  as.matrix(rownames.value = snplist)


# run MR
protein <- snakemake@wildcards[["protein"]]

df_res <- data_MR[Protein == protein] %>% 
  group_by(r2_thresh, P_thresh, Protein) %>% 
  nest() %>% 
  # run MR per instrument set
  mutate(MR = map(data, run_MR_all, ldrho)) %>% 
  select(-data) %>% 
  unnest(cols = MR)

fwrite(df_res, snakemake@output[[1]])

