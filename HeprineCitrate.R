library(xlsx)
setwd("/home/rstudio/host/beef_git_tensor/Beef-Metabolomics/HeprineCitrate.xlsx")
file1 <- "/home/rstudio/host/beef_git_tensor/Beef-Metabolomics/HeprineCitrate.xlsx"
file2 <- "/home/rstudio/host/beef_git_tensor/Beef-Metabolomics/HeprineCitrate2.xls"

key_string1 = readxl::read_excel(path = file1, sheet = "HeprineCitrate_KS")
key_string2 = readxl::read_excel(path = file2, sheet = "1", col_names = TRUE)
key_string3 = readxl::read_excel(path = file2, sheet = "2", col_names = TRUE)

db = rbind(key_string2,key_string3)


final = left_join(key_string1,db)

colnames(final)

openxlsx::write.xlsx(final,"Heprine_Citrate_Sub.xlsx" )