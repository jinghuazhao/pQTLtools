source("workflow/r/MR_functions.R")

protein <- snakemake@wildcards[["protein"]]
snplist <- read_lines(snakemake@input[["snplist"]])
ld <- fread(snakemake@input[["ld"]], col.names = snplist) %>%
      as.matrix(rownames.value = snplist)

df_res <- fread(snakemake@input[["data_MR"]]) %>%
          setnames(c("rsID", "beta_Prot", "se_Prot", "beta_HF", "se_HF"),
                   c("SNP",  "beta.x",    "se.x",    "beta.y",  "se.y")) %>%
          filter(Protein == protein) %>%
          group_by(r2_thresh, P_thresh, Protein) %>% 
          nest() %>% 
          mutate(MR = map(data, run_MR, ld)) %>%
          select(-data) %>% 
          unnest(cols = MR)
fwrite(df_res, snakemake@output[[1]])
