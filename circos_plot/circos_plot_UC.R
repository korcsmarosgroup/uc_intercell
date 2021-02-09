#Visualizing condition specific ligand-receptor connection between myofibroblast and regulator T cells

#Install circlize package if needed
#install.packages("circlize")

#call package
library(circlize)

#setting up working directory
setwd('/Users/lgul/Documents/OneDrive_Norwich_BioScience_Institutes/sc_data_analysis/OmniPath2/minor_revision/UC')

#!/usr/bin/env Rscript
require(circlize)
mat <- read.csv('myofibroblast_Treg_circos_table_v5.csv')
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,-1])
png('myo-treg_v3.png', width = 7, height = 7, units = 'in', res = 300)

#coloring edges based on pathways of receptor
grid.col = c(ITGA = '#FF6C69FF',CX3CR1 = '#FF6C69FF', PTPR = '#FF6C69FF', CCR = '#FF6C69FF', SIGLEC10 = '#FF6C69FF',TREM2 = '#FF6C69FF', SELP = '#FF6C69FF', FPR1 = '#FF6C69FF', KIR3DL2 = '#FF6C69FF', MADCAM1 = '#FF6C69FF',
             KLRC1 = '#EAD800FF', IL1B = '#EAD800FF', 
             GPR160 = 'grey', F3 = 'grey', MCAM = 'grey', FAT4 = 'grey', XPR1 = 'grey', RAMP2 = 'grey', CNR2 = 'grey', 
             DDR1 = '#D47A02FF', PDGFRA = '#D47A02FF', EPHA2 = '#D47A02FF', INSR = '#D47A02FF', FGFR3 = '#D47A02FF', PTPRF = '#D47A02FF', VANGL1 = '#D47A02FF', IGF1R = '#D47A02FF', NRP1 = '#D47A02FF', S1PR1 = '#D47A02FF', 
             ACVRL1 = '#B3000DFF', 
             TLR2 = '#FD9409FF', CXCL9 = '#FD9409FF', TNFRSF8 = '#FD9409FF', LGR4 = '#FD9409FF',
             LRP6 = '#008A5BFF', PTK7 = '#008A5BFF',
             DLK2 = 'blue',
             ADM = 'black',
             ANGPT1 = 'black',
             ANPEP = 'black',
             ANXA1 = 'black',
             BDNF = 'black',
             BMP2 = 'black',
             CAV1 = 'black',
             CCLs = 'black',
             CD14 = 'black',
             CD52 = 'black',
             COL21A1 = 'black',
             COL6A3 = 'black',
             CSF1 = 'black',
             CTGF = 'black',
             CX3CL1 = 'black',
             CXCLs = 'black',
             DCHS1 = 'black',
             DDR1 = 'black',
             DEFB1 = 'black',
             EDN1 = 'black',
             EFNAs = 'black',
             EPHA4 = 'black',
             FGFs = 'black',
             FLT1 = 'black',
             FLT3LG = 'black',
             GNAI2 = 'black',
             GPC3 = 'black',
             HLAs = 'black',
             HRAS = 'black',
             IAPP = 'black',
             IL10 = 'black',
             IL12RB1 = 'black',
             IL1R1 = 'black',
             ITGAV = 'black',
             ITGB3 = 'black',
             KITLG = 'black',
             LAMAs = 'black',
             LAMBs = 'black',
             LAMC1 = 'black',
             LGALS3BP = 'black',
             LRP1 = 'black',
             MCAM = 'black',
             NID1 = 'black',
             NOTCH1 = 'black',
             NRG1 = 'black',
             NRP1 = 'black',
             NTF3 = 'black',
             PDGFs = 'black',
             PDGFRA = 'black',
             PGF = 'black',
             POMC = 'black',
             PTHLH = 'black',
             PTN = 'black',
             PTPRJ = 'black',
             PTPRM = 'black',
             RSPOs = 'black',
             SCT = 'black',
             SDC2 = 'black',
             SEMA3s = 'black',
             TGFBs = 'black',
             THBS1 = 'black',
             TNC = 'black',
             TNF = 'black',
             TNFSF12 = 'black',
             VCAN = 'black',
             VEGFs= 'black',
             VWF= 'black',
             WNTs = 'black')

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