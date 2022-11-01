# -----------------index data analysis-----------------
# ---------------------SimchaUbo-----------------------
# ---------------------2022.10.18----------------------



# transform index.sorting.fcs file to MFI data matrix---------------------------------->> 
rm(list = ls())
gc()
# BiocManager::install('flowCore')
# devtools::install_github("Kawameicha/indexSort")

library(tidyverse)
library(flowCore)
library(indexSort)

inputFCS.lung <- read.FCS("E:/Results/FCM/Mouse/Mouse_project/Itga4/20220816_W4_Sort/20220816_W4_SORT_INX_Lung3_001_018.fcs")
inputFCS.liver <- read.FCS("E:/Results/FCM/Mouse/Mouse_project/Itga4/20220816_W4_Sort/20220816_W4_SORT_INX_Liver3_001_023.fcs")
inputFCS.kidney <- read.FCS("E:/Results/FCM/Mouse/Mouse_project/Itga4/20220816_W4_Sort/20220816_W4_SORT_INX_Kidney2_001_026.fcs")


result.lung <- retrieve_index(inputFCS.lung)
result.liver <- retrieve_index(inputFCS.liver)
result.kidney <- retrieve_index(inputFCS.kidney)
# lung_d3_Plate1_FACS data
write.table(result.lung,'E:/Results/FCM/Mouse/Mouse_project/Itga4/lung_d3_1.txt',quote=F, sep = '\t')

# liver_d3_Plate1_FACS data
write.table(result.liver,'E:/Results/FCM/Mouse/Mouse_project/Itga4/liver_d3_1.txt',quote=F, sep = '\t')

# kidney_d2_Plate1_FACS data
write.table(result.kidney,'E:/Results/FCM/Mouse/Mouse_project/Itga4/kidney_d2_1.txt',quote=F, sep = '\t')

explore_plate(result.kidney, param = c("APC.A", "PE.A", "DAPI.A", "PE.Cy7.A", "PerCP.Cy5.5.A", 'APC.Cy7.A','Time'), dim = c(16,24),circle = 6)

# read data----------------------------------------------------------------------------->> 
source('/work/mnt/sourceR/Seurat_utils.R')
setwd('/work/mnt/00mouse_data/v8/T_submc/')

seu1 = readRDS('saved_work/submc_T_v2.rds') # object with annot 
seu2 = readRDS('output/merge_T.rds') # object without annot

ab3369 = read.table('saved_work/SB50/AB3369.txt') 
ab3370 = read.table('saved_work/SB50/AB3370.txt') 
ab3371 = read.table('saved_work/SB50/AB3371.txt') 

# add annot to merge data--------------------------------------------------------------->>  
seu2@meta.data[,'cell_id'] = rownames(seu2@meta.data)

tmp1 = seu1@meta.data %>% dplyr::select(cell_id,annots)
tmp2 = seu2@meta.data
tmp = merge(tmp2,tmp1,by = 'cell_id')
seu2@meta.data$annots = ''
seu2@meta.data[intersect(rownames(seu1@meta.data),rownames(seu2@meta.data)),'annots'] = seu1@meta.data[intersect(rownames(seu1@meta.data),rownames(seu2@meta.data)),'annots']


# extract sb50 data to compare project location vs rest data---------------------------->>  


sb50 <- subset(seu2, cell_id %in% cell_show)
sb50_cut = subset(sb50,nCount_RNA >=450)

sb50 <- NormalizeData(sb50, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
sb50_cut <- NormalizeData(sb50_cut, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)


# form well-coordinate--------------------------------------------------------------->>  
# create a MARS-dependent well index
M1 = matrix(seq(1,384,4),8,12,byrow = TRUE)
M2 = matrix(seq(2,384,4),8,12,byrow = TRUE)
M3 = matrix(seq(3,384,4),8,12,byrow = TRUE)
M4 = matrix(seq(4,384,4),8,12,byrow = TRUE)

colindex = seq(1:24)
colindex[c(TRUE,FALSE)] = seq(1:12)
colindex[c(FALSE,TRUE)] = seq(13,24,1)

rowindex = seq(1:16)
rowindex[c(TRUE,FALSE)] = seq(1:8)
rowindex[c(FALSE,TRUE)] = seq(9,16,1)


M12 = cbind(M1,M2)
M12 = M12[,colindex]

M34 = cbind(M3,M4)
M34 = M34[,colindex]

M = rbind(M12,M34)

M = M[rowindex,]

M_vec = as.vector(t(M))

M_vec

# standard 384 index
F_vec = sapply(LETTERS[1:16], function(x) paste0(x,1:24))
F_vec = as.vector(F_vec)
F_vec   

map = as.data.frame(cbind(index = F_vec, well = paste0('W',M_vec+1)))

map$row_index = str_split(map$index,pattern = '',simplify = TRUE,n=2)[,1] %>% as.integer()
map$col_index = str_split(map$index,pattern = '',simplify = TRUE,n=2)[,2] %>% as.integer()
map$well_index = str_split(map$well,pattern = '', simplify = TRUE, n = 2)[,2] %>% as.integer()

map <- map%>% arrange(well_index)


# arrange fj_data--------------------------------------------------------------->>  


fj_ab3369$Idx = paste0(fj_ab3369$IdxRow,fj_ab3369$IdxCol)
rownames(fj_ab3369) <- NULL
fj_ab3369 = fj_ab3369[match(map$index, fj_ab3369$Idx),]
removeRowsAllNa(fj_ab3369) -> fj_ab3369
rownames(fj_ab3369) = paste0('AB3369.W',c(2:15,18:31,34:385)) # rm O1,O2,P1,P2 index
fj_ab3369$tissue = 'liver'

fj_ab3370$Idx = paste0(fj_ab3370$IdxRow,fj_ab3370$IdxCol)
rownames(fj_ab3370) <- NULL
fj_ab3370 = fj_ab3370[match(map$index, fj_ab3370$Idx),]
removeRowsAllNa(fj_ab3370) -> fj_ab3370
rownames(fj_ab3370) = paste0('AB3370.W',c(386:399,402:415,418:769)) # rm O1,O2,P1,P2 index
fj_ab3370$tissue = 'lung'

fj_ab3371$Idx = paste0(fj_ab3371$IdxRow,fj_ab3371$IdxCol)
rownames(fj_ab3371) <- NULL
fj_ab3371 = fj_ab3371[match(map$index, fj_ab3371$Idx),]
removeRowsAllNa(fj_ab3371) -> fj_ab3371
rownames(fj_ab3371) = apaste0('AB3371.W',c(770:783,786:799,802:886,888:1153)) # lack of A8 index data
fj_ab3371$tissue = 'kidney'

fj_data = rbind(fj_ab3369,fj_ab3370,fj_ab3371)
#dim(fj_data)
head(fj_data)


# plot umi-normalization and FACS MFI----------------------------------------------->>  
# base plot
plot_FACS_scRNA <- function(sc_object, fj_object, gene_id, flow_id, title = '')
    { 
    plot(sc_object@assays$RNA@data[gene_id,],fj_object[rownames(sc_object@meta.data),flow_id],main = title,xlab = '', ylab='')
    }

# ggplot
plot_FACS_scRNA2 <- function(sc_object, fj_object, gene_id, flow_id, title = '', na.omit = FALSE)
    { 
    df = as.data.frame(
                        cbind(
                        umi_norm = sc_object@assays$RNA@data[gene_id,], 
                        MFI = fj_object[rownames(sc_object@meta.data),flow_id],
                        tissue = fj_object[rownames(sc_object@meta.data),'tissue']        
                             )
                      )
    df$umi_norm = as.numeric(df$umi_norm)
    df$MFI = as.numeric(df$MFI)

    if(na.omit) {df = subset(df,!is.na(MFI)); df = subset(df, umi_norm!=0)}
    
    ggplot(df,aes(umi_norm,MFI,color = tissue))+
        geom_point()+
        geom_smooth(method = "lm")+
        ggtitle(title)
    
    }


options(repr.plot.height = 6, repr.plot.width = 8)
p0 <- plot_FACS_scRNA2(sb50_cut, fj_data,'Cd8a','PE.Cy7.A',title='CD8-Cd8a',na.omit = TRUE)

p1 <- plot_FACS_scRNA2(sb50_cut, fj_data,'Itga4','APC.A',title='CD49d-Itga4',na.omit = TRUE)

p2 <- plot_FACS_scRNA2(sb50_cut, fj_data,'Klrg1','PE.A',title='KLRG1-Klrg1',na.omit = TRUE)

p3 <- plot_FACS_scRNA2(sb50_cut, fj_data,'Ptprc','APC.Cy7.A',title='CD45-Ptprc',na.omit = TRUE)

(p1 +p2) /(p3 +p0)	
