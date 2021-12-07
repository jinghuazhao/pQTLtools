# SomaScan V4.1

HOME <- Sys.getenv("HOME")
SomaScanV4.1 <- openxlsx::read.xlsx(file.path(HOME,"SomaLogic","doc","v4.1","SomaScan v4.1 Protein Content Menu.xlsx"),
                                    sheet="SomaScanV4.1_Protein_Content", colNames=TRUE, skipEmptyRows=TRUE, cols=1:6, startRow=3)
save(SomaScanV4.1,file='SomaScanV4.1.rda',compress='xz')
