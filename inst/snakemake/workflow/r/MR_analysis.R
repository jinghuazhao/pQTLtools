source("workflow/r/MR_functions.R")

protein <- snakemake@wildcards[["protein"]]
trait <- snakemake@wildcards[["trait"]]
snplist <- read_lines(snakemake@input[["snplist"]])
ld <- fread(snakemake@input[["ld"]], col.names = snplist) %>%
      as.matrix(rownames.value = snplist)

df_res <- fread(snakemake@input[["data_MR"]]) %>%
          rename("SNP"="rsID", "beta.x"="beta_Prot", "se.x"="se_Prot", "beta.y"=paste0("beta_",trait),"se.y"= paste0("se_",trait)) %>%
          filter(Protein == protein) %>%
          group_by(r2_thresh, P_thresh, Protein) %>% 
          nest() %>% 
          mutate(MR = map(data, run_MR, ld)) %>%
          select(-data) %>% 
          unnest(cols = MR)
fwrite(df_res, snakemake@output[[1]])
