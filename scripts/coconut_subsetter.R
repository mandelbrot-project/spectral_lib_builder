library(dplyr)
library(readr)


np_atlas_2019_12 <- read_delim("is_fragmentation/npatlas_data/np_atlas_2019_12.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

npatlas_for_frag <- np_atlas_2019_12 %>% 
            select(NPAID, SMILES)


coconut_for_frag <- COCONUT4MetFrag %>% 
                select(Identifier, SMILES)



write.table(coconut_for_frag, file='/Volumes/COMMON FASIE-FATHO/Coconut_Frag/coconut_for_frag.txt', row.names=FALSE, quote=FALSE, sep=" ")
write.table(npatlas_for_frag, file='~/is_fragmentation/npatlas_data/npatlas_for_frag.txt', row.names=FALSE, quote=FALSE, sep=" ")
