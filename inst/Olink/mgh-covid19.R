
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(OlinkAnalyze))

document_example <- function()
# https://github.com/Olink-Proteomics/OlinkRPackage
{
  # visualize the NPX distribution per sample per panel, example for one panel
  olink_dist_plot(npx_data1 %>% filter(Panel == 'Olink Cardiometabolic')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c('turquoise3', 'red'))

  olink_qc_plot(npx_data1 %>% filter(Panel == 'Olink Cardiometabolic')) +
  scale_color_manual(values = c('turquoise3', 'red'))

  # identify bridge samples
  bridge_samples <- intersect(x = npx_data1$SampleID,
                              y = npx_data2$SampleID)

  # bridge normalize
  bridge_normalized_data <- olink_normalization(df1 = npx_data1,
                            df2 = npx_data2,
                            overlapping_samples_df1 = bridge_samples,
                            df1_project_nr = "20200001",
                            df2_project_nr = "20200002",
                            reference_project = "20200001")
  # t-test npx_data1
  ttest_results_NPX1 <- olink_ttest(df = npx_data1,
                                    variable = "Treatment")

  # select names of the top #10 most significant proteins
  ttest_sign_NPX1 <- ttest_results_NPX1 %>%
  head(n=10) %>%
  pull(OlinkID)

  # volcano plot with annotated top #10 most significant proteins
  olink_volcano_plot(p.val_tbl = ttest_results_NPX1,
                     olinkid_list = ttest_sign_NPX1) +
  scale_color_manual(values = c('turquoise3', 'red'))
}

# MGH data

# https://www.olink.com/application/mgh-covid-19-study/

HOME <- Sys.getenv("HOME")
dir <- file.path(HOME,"R","work","MGH","MGH_Olink_COVID_Apr_27_2021/")

NPX_data <- read_NPX(filename = file.path(dir,"MGH_COVID_OLINK_NPX.txt"))

# visualize the NPX distribution per sample per panel, example for one panel
olink_dist_plot(NPX_data %>% filter(Panel == 'Cardiometabolic')) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c('turquoise3', 'red'))

olink_qc_plot(NPX_data %>% filter(Panel == 'Cardiometabolic')) +
  scale_color_manual(values = c('turquoise3', 'red'))
