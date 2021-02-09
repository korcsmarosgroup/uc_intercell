#Visualizing condition specific ligand-receptor connection between myofibroblast and regulator T cells
#Authors: Lejla Gul, Denes Turei

#Input: Table describing the interactions
#Output: PDF file

#Install circlize package if needed
#install.packages("circlize")

#call the package
library(circlize)

#settin up working directory
setwd('/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/OmniPath2/minor_revision')

#!/usr/bin/env Rscript
require(circlize)

#Import input file
mat <- read.csv('myofibroblast_Treg_circos_table_H_v2.csv')
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,-1])
pdf('myo-treg_H.png', width = 7, height = 7, units = 'in', res = 300)

#coloring edges based on pathways of receptor
#black = ligand
#grey = 'other' pathway
# #B3000DFF = TGF-beta pathway
# #3DC041FF = Hedgehog pathway
# #EAD800FF = JAK/STAT pathway
# #D47A02FF = RTK pathway
# #FF6C69FF = Innate immune pathways
# #FD9409FF = TLR pathway
# #008A5BFF = WNT pathway

grid.col = c(ADM = 'black', AGRN = 'black',
             ANGPT1 = 'black',
             ANXA1 = 'black',
             APOs = 'black',
             APP = 'black',
             BDNF = 'black',
             BMPs = 'black',
             BMPR2 = 'black',
             CALR = 'black',
             CCLs = 'black',
             CD14 = 'black',
             CD47 = 'black',
             CD81 = 'black',
             COL6A3 = 'black',
             COPA = 'black',
             CSF1 = 'black',
             CTGF = 'black',
             CX3CL1 = 'black',
             CXCLs = 'black',
             DCBLD2 = 'black',
             ECM1 = 'black',
             EDN1 = 'black',
             EFNAs = 'black',
             EPHB4 = 'black',
             FGFs = 'black',
             FLT3LG = 'black',
             GDF11 = 'black',
             GDNF = 'black',
             GNAI2 = 'black',
             GPC3 = 'black',
             HRAS = 'black',
             IAPP = 'black',
             ICAMs = 'black',
             ILs = 'black',
             ITGA4 = 'black',
             ITGB1 = 'black',
             KITLG = 'black',
             LAMAs = 'black',
             LAMBs = 'black',
             LAMC1 = 'black',
             LGALS9 = 'black',
             MDK = 'black',
             NCSTN = 'black',
             NOTCH1 = 'black',
             NRG1 = 'black',
             NRP2 = 'black',
             NTF3 = 'black',
             PDGFs = 'black',
             PDGFRA = 'black',
             PGF = 'black',
             PLXNB2 = 'black',
             POMC = 'black',
             PTHLH = 'black',
             PTN = 'black',
             PTPRJ = 'black',
             RAC1 = 'black',
             RPSA = 'black',
             RTN4 = 'black',
             SCT = 'black',
             SEMA3s = 'black',
             TGFB1 = 'black',
             THBS1 = 'black',
             TNC = 'black',
             TNF = 'black',
             TNFSF12 = 'black',
             VCAN = 'black',
             VEGFs = 'black',
             VWF = 'black',
             WNTs = 'black',
             ADCY9 = 'grey',
             CD109 = '#B3000DFF',
             GPC1 = 'grey',
             LRP = 'grey',
             P2RY6 = 'grey',
             PLXNA3 = 'grey',
             RAMP3 = 'grey',
             SEMA4 = 'grey',
             GPC = '#3DC041FF',
             ITGA = '#FF6C69FF',
             ACKR3 = '#FF6C69FF',
             CCR = '#FF6C69FF',
             CX3CL1 = '#FF6C69FF',
             GP1BA = '#FF6C69FF',
             PTPR = '#FF6C69FF',
             CSF2 = '#EAD800FF',
             FLT4 = '#D47A02FF',
             NGFR = '#D47A02FF',
             SCARF1 = '#FF6C69FF',
             THBS1 = '#FF6C69FF',
             LTK = '#D47A02FF',
             PDGFB = '#D47A02FF',
             SDC = '#D47A02FF',
             BMPR1A = '#B3000DFF',
             ITGB8 = '#FF6C69FF',
             VCAM1 = '#FF6C69FF',
             CXCL10 = '#FD9409FF',
             IL1RN = '#EAD800FF',
             TLR = '#FD9409FF',
             TSHR = '#FD9409FF',
             CD36 = '#FF6C69FF',
             FZD = '#008A5BFF')

chordDiagram(mat, grid.col = grid.col, annotationTrack = 'grid', preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  ylim = c(0, .3),
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      mean(xlim),
      -ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, .5),
      cex = .6
    )
  },
  bg.border = NA
)
circos.clear()
dev.off()