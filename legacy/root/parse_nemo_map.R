library(data.table)
map_nemo <- fread("./Nemo3/NOBGS 3/r-e5/GENO/test-Chr1_V1_c1_env1_1_1.map")
GTs <- fread("./Nemo3/NOBGS 3/r-e5/GENO/test-Chr1_V1_c1_env1_1_1.snp_geno")

sample_info <- GTs[,.(pop,ID)]
GTs <- as.matrix(GTs[1:.N,6:ncol(GTs),with=FALSE])


# original recombination map
map <- readRDS("./Nemo3/map_100kb_100qtn_100delet/chromosome_maps_100kb_rds/ rec_map 1 .rds")
map[,indx := .I]

# parse nemo map
nemo_map <- data.table(marker=map_nemo$trait.locus,do.call(rbind,strsplit(map_nemo$trait.locus,".",fixed=TRUE)))
setnames(nemo_map, c("V1","V2"),c("type","idx"))


nemo_map[,idx:=as.numeric(idx)+1] ## because nemo starts from 0

# indexes in Nemo map
ntrl_idx   <- nemo_map[type=="ntrl",idx]
quanti_idx <- nemo_map[type=="quant",idx]
delet_idx <- nemo_map[type=="delet" ,idx] 

# indexes in original map
indx_ntrl   <- map[type=="ntrl"][ntrl_idx,indx]
indx_quanti <- map[type=="QTN"][quanti_idx,indx]
indx_delet <- map[type=="delet"][delet_idx,indx]

# join
map[indx_ntrl,nemo_marker     := nemo_map[type=="ntrl",marker]]
map[indx_quanti,nemo_marker   := nemo_map[type=="quant",marker]]
map[indx_delet,nemo_marker   := nemo_map[type=="delet" ,marker]]

# match and reorder
map <- map[nemo_marker %in% colnames(GTs)]
GTs <- GTs[,map$nemo_marker]
map[,Pos:=bp]
map[,marker := paste(paste0("Chr",Chr),Pos,sep=":")]
colnames(GTs) <- map$marker
setorder(map,Chr,bp)

GTs <- GTs[,map$marker]
map <- map[,.(Chr,Pos,marker,type)]
