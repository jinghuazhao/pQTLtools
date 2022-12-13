suppressMessages(library(tidyverse))
suppressMessages(library(metafor))

df_obs <- read_csv(snakemake@input[[1]])
run_meta <- function(df, y_col = beta, se_col = se, method = "FE", ...,
                     report=c("beta", "se", "pval", "QEp"))
{
  y_col <- enquo(y_col)
  se_col <- enquo(se_col)
  betas <- pull(df, !!y_col)
  ses <- pull(df, !!se_col)
  meta_res <- tryCatch(rma(yi = betas, sei = ses, method = method), error = function(e) NULL)
  df.res <- if (!is.null(meta_res)){
    map(report, ~`[[`(meta_res,.) %>% as.numeric()) %>%
    set_names(report) %>%
    as_tibble()
  } else rep(NaN, length(report)) %>% as.list %>% set_names(report) %>% as_tibble()
}

df_meta <- df_obs %>% 
  group_by(Protein) %>% 
  nest() %>% 
  mutate(meta = map(data, run_meta)) %>% 
  select(-data) %>% 
  unnest() %>% 
  mutate(RR = exp(beta),
         RR_LCI =  exp(beta - qnorm(0.975)*se),
         RR_UCI = exp(beta + qnorm(0.975)*se))
write_csv(df_meta, snakemake@output[[1]])
