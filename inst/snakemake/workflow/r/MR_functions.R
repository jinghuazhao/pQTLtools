# Helper function to run MR

library(data.table)
library(tidyverse)
library(MendelianRandomization)

# Install MR package
MR_package_ver <- snakemake@config[["MR_package_ver"]]
install_MR <- function(ver){
  if (MR_package_ver == "latest")
    install.packages("MendelianRandomization")
  else {
    URL <- paste0("https://cran.r-project.org/src/contrib/Archive/MendelianRandomization/MendelianRandomization_",
                MR_package_ver, ".tar.gz")
    install.packages(URL, repos = NULL, type = "source")
  }
}

if (!require("MendelianRandomization")){
  install_MR()
} else if(packageVersion("MendelianRandomization") != MR_package_ver) {
  install_MR()
}


# Blank results
error_df <- function(df_mr){
  n_ins <- nrow(df_mr)
  data.table(N_ins = n_ins,
             beta = NaN, SE = NaN,
             LCI = NaN, UCI = NaN,
             P_value = NaN)
}

blank_df <- function(){
  data.table(N_ins = 0,
             beta = NA, SE = NA,
             LCI = NA, UCI = NA,
             P_value = NA, Method = NA)
}

# Single instrument MR (ratio method)
single.ins.MR <- function (df_mr, harmonise=F){
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
  return(res)
}

# MR IVW using MendelianRandomization package
# (more flexible as additional parameters are allowed, e.g. random-effect model)
# MR IVW w/ correlation using MendelianRandomization package
MR_IVW <- function (df_mr, ldrho=matrix(), harmonise=F,  ...) {
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]

    # add correlation matrix to the model if LD rho is supplied
    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ldrho,
                         snps = df_mr[,SNP])
    MR.output <- mr_ivw(MR.input, correl = TRUE, ...)

    data.table(N_ins = MR.output@SNPs,
               beta = MR.output@Estimate,
               SE = MR.output@StdError,
               LCI = MR.output@CILower,
               UCI = MR.output@CIUpper,
               P_value = MR.output@Pvalue)
  },
  error = function(e) error_df(df_mr)
  )
  return(res)
}


# MR Egger
MR_Egger <- function (df_mr, ldrho=matrix(), harmonise=F, ...) {
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]

    MR.input <- mr_input(bx = df_mr[,beta.x],
                         bxse = df_mr[,se.x],
                         by = df_mr[,beta.y],
                         byse = df_mr[,se.y],
                         corr = ldrho,
                         snps = df_mr[,SNP])

    MR.output <- mr_egger(MR.input, correl = TRUE, ...)

    data.table(N_ins = MR.output@SNPs,
               beta = c(MR.output@Estimate, MR.output@Intercept),
               SE = c(MR.output@StdError.Est, MR.output@StdError.Int),
               LCI = c(MR.output@CILower.Est, MR.output@CILower.Int),
               UCI = c(MR.output@CIUpper.Est, MR.output@CIUpper.Int),
               P_value = c(MR.output@Pvalue.Est, MR.output@Pvalue.Int)
    )
  },
  error = function(e) rbind(error_df(df_mr), error_df(df_mr))
  )

  res[, Method := c("MR_Egger", "EggerInt")]
  return(res)
}


# MR IVW-PCA method
# var_exp = expected variance in the risk factor explained by principal components
MR_PCA <- function(df_mr, ldrho, harmonise=F, var_exp=0.99){
  res <- tryCatch({
    if (harmonise) df_mr[, beta.y := ifelse(A1.x==A1.y, beta.y, -beta.y)]

    attach(df_mr)
    Phi = (beta.x / se.y) %o% (beta.x / se.y) * ldrho
    # summary(prcomp(Phi, scale=FALSE))

    K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2 /
                       sum((prcomp(Phi, scale=FALSE)$sdev^2)))
              > var_exp)[1]
    # K is number of principal components to include in analysis
    # this code includes principal components to explain 99% of variance in the risk factor

    betaXG0 = as.numeric(beta.x%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
    betaYG0 = as.numeric(beta.y%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])

    Omega = se.y %o% se.y * ldrho

    pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]

    #Calculate the MR beta and se
    beta_IVWcorrel.pc <- solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
    beta_IVWcorrel.pc <- as.numeric(beta_IVWcorrel.pc)
    se_IVWcorrel.fixed.pc <- sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)) %>%
      as.numeric
    #Z-score and P-value
    Z_score = beta_IVWcorrel.pc / se_IVWcorrel.fixed.pc
    P_value = 2*pnorm(-abs(Z_score))
    LCI <- beta_IVWcorrel.pc - qnorm(0.975)*se_IVWcorrel.fixed.pc
    UCI <- beta_IVWcorrel.pc + qnorm(0.975)*se_IVWcorrel.fixed.pc

    detach(df_mr)

    data.table(N_ins = nrow(df_mr),
               beta = beta_IVWcorrel.pc,
               SE = se_IVWcorrel.fixed.pc,
               LCI = LCI, UCI = UCI,
               P_value = P_value)
  },
  error = function(e) error_df(df_mr)
  )

  return (res)
}

# Wrapper to run all MR based on n_instrument
run_MR_all <- function(df_mr, ldrho){
  setDT(df_mr)
  if (nrow(df_mr) == 0) {
    res <- blank_df()
  } else if (nrow(df_mr) == 1) {
    res <- single.ins.MR(df_mr)
    res[, `:=`(Method = "Wald_1Ins")]
  } else {

    insSNPs <- df_mr$SNP

    ldrho_ins <- ldrho[insSNPs, insSNPs]

    # Exclude SNP if LDcorr results in NaN
    NaN_ins <- which(is.nan(ldrho_ins), T) %>% rownames
    if (!is.null(NaN_ins)){
      ins <- rownames(ldrho_ins) %>% .[which(!. %in% NaN_ins)]
      ldrho_ins <- ldrho_ins[ins, ins]

      df_mr <- df_mr[SNP %in% ins]
    }

    res_IVW <- MR_IVW(df_mr, ldrho=ldrho_ins, harmonise = F)
    res_IVW[, Method := "IVW"]
    res_Egger <- MR_Egger(df_mr, ldrho=ldrho_ins, harmonise = F)

    res_PCA_0.99 <- MR_PCA(df_mr, ldrho_ins, harmonise = F, var_exp=0.99) %>%
      .[, `:=`(Method = "PCA_0.99")]

    res_PCA_0.90 <- MR_PCA(df_mr, ldrho_ins, harmonise = F, var_exp=0.90) %>%
      .[, `:=`(Method = "PCA_0.90")]

    res <- rbind(res_IVW, res_PCA_0.90, res_PCA_0.99, res_Egger)
  }
  return(res)
}
