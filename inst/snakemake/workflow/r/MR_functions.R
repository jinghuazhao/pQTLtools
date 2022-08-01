suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(MendelianRandomization))

MR_package_ver <- snakemake@config[["MR_package_ver"]]

install_MR <- function()
{
  if (MR_package_ver == "latest") install.packages("MendelianRandomization")
  else {
    URL <- paste0("https://cran.r-project.org/src/contrib/Archive/MendelianRandomization/MendelianRandomization_", MR_package_ver, ".tar.gz")
    install.packages(URL, repos = NULL, type = "source")
  }
}

if (!require("MendelianRandomization") | packageVersion("MendelianRandomization") != MR_package_ver) install_MR()

error_df <- function(df_mr)
  data.table(N_ins = nrow(df_mr),
             beta = NaN, SE = NaN,
             LCI = NaN, UCI = NaN,
             P_value = NaN)

blank_df <- function()
  data.table(N_ins = 0,
             beta = NA, SE = NA,
             LCI = NA, UCI = NA,
             P_value = NA, Method = NA)

single.ins.MR <- function (df_mr, harmonise=FALSE)
{
  if (nrow(df_mr) != 1) stop ("N instrument is not 1")
  if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]
  df_mr[, `:=`(beta.mr = beta.y / beta.x,
               se.mr = sqrt((se.y^2/beta.x^2) + (beta.y^2 *se.x^2/beta.x^4))
  )][,`:=`(lci.mr = beta.mr - qnorm(0.975)*se.mr,
           uci.mr = beta.mr + qnorm(0.975)*se.mr,
           P.mr = 2*pnorm(-abs(beta.mr / se.mr))
  )]
  res <- data.table(N_ins = 1,
                    beta = df_mr[,beta.mr],
                    SE = df_mr[,se.mr],
                    LCI = df_mr[,lci.mr],
                    UCI = df_mr[,uci.mr],
                    P_value = df_mr[,P.mr])
}

MR_IVW <- function (df_mr, ld=matrix(), harmonise=FALSE,  ...)
{
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]
    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ld,
                         snps = df_mr[,SNP])
    MR.output <- mr_ivw(MR.input, correl = TRUE, ...)
    data.table(N_ins = MR.output@SNPs,
               beta = MR.output@Estimate,
               SE = MR.output@StdError,
               LCI = MR.output@CILower,
               UCI = MR.output@CIUpper,
               P_value = MR.output@Pvalue)
  }, error = function(e) error_df(df_mr))
}

MR_Egger <- function (df_mr, ld=matrix(), harmonise=FALSE, ...)
{
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]
    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ld,
                         snps = df_mr[,SNP])
    MR.output <- mr_egger(MR.input, correl = TRUE, ...)
    data.table(N_ins = MR.output@SNPs,
               beta = c(MR.output@Estimate, MR.output@Intercept),
               SE = c(MR.output@StdError.Est, MR.output@StdError.Int),
               LCI = c(MR.output@CILower.Est, MR.output@CILower.Int),
               UCI = c(MR.output@CIUpper.Est, MR.output@CIUpper.Int),
               P_value = c(MR.output@Pvalue.Est, MR.output@Pvalue.Int)
    )
  }, error = function(e) rbind(error_df(df_mr), error_df(df_mr)))
  res[, Method := c("MR_Egger", "EggerInt")]
  res
}

MR_PCA <- function(df_mr, ld, harmonise=FALSE, var_exp=0.99)
{
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]
    attach(df_mr)
    Phi <- (beta.x / se.y) %o% (beta.x / se.y) * ld
    K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2 / sum((prcomp(Phi, scale=FALSE)$sdev^2))) > var_exp)[1]
    betaXG0 <- as.numeric(beta.x%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
    betaYG0 <- as.numeric(beta.y%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
    Omega <- se.y %o% se.y * ld
    pcOmega <- t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]
    beta_IVWcorrel.pc <- solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
    beta_IVWcorrel.pc <- as.numeric(beta_IVWcorrel.pc)
    se_IVWcorrel.fixed.pc <- sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)) %>% as.numeric()
    Z_score <- beta_IVWcorrel.pc / se_IVWcorrel.fixed.pc
    P_value <- 2*pnorm(-abs(Z_score))
    LCI <- beta_IVWcorrel.pc - qnorm(0.975)*se_IVWcorrel.fixed.pc
    UCI <- beta_IVWcorrel.pc + qnorm(0.975)*se_IVWcorrel.fixed.pc
    detach(df_mr)
    data.table(N_ins = nrow(df_mr),
               beta = beta_IVWcorrel.pc,
               SE = se_IVWcorrel.fixed.pc,
               LCI = LCI, UCI = UCI,
               P_value = P_value)
  }, error = function(e) error_df(df_mr))
}

run_MR_all <- function(df_mr, ld)
{
  setDT(df_mr)
  if (nrow(df_mr) == 0) res <- blank_df()
  else if (nrow(df_mr) == 1) {
    res <- single.ins.MR(df_mr)
    res[, `:=`(Method = "Wald_1Ins")]
  } else {
    insSNPs <- df_mr$SNP
    ld_ins <- ld[insSNPs, insSNPs]
    NaN_ins <- which(is.nan(ld_ins), T) %>% rownames
    if (!is.null(NaN_ins)){
      ins <- rownames(ld_ins) %>% .[which(!. %in% NaN_ins)]
      ld_ins <- ld_ins[ins, ins]
      df_mr <- df_mr[SNP %in% ins]
    }
    res_IVW <- MR_IVW(df_mr, ld=ld_ins, harmonise = FALSE)
    res_IVW[, Method := "IVW"]
    res_Egger <- MR_Egger(df_mr, ld=ld_ins, harmonise = FALSE)
    res_PCA_0.99 <- MR_PCA(df_mr, ld_ins, harmonise = FALSE, var_exp=0.99) %>%
      .[, `:=`(Method = "PCA_0.99")]
    res_PCA_0.90 <- MR_PCA(df_mr, ld_ins, harmonise = FALSE, var_exp=0.90) %>%
      .[, `:=`(Method = "PCA_0.90")]
    res <- rbind(res_IVW, res_PCA_0.90, res_PCA_0.99, res_Egger)
  }
  res
}
