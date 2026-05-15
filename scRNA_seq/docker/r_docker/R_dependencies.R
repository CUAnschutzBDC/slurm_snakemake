# Note Rhtslib is often an issue, renv::install("bioc::Rhtslib", type = "binary")
# most recently worked
library(tidyverse)
library(BiocManager)
library(Seurat)
library(biomaRt)
library(scAnalysisR) # CUAnschutzBDC/scAnalysisR
library(clustifyr)
library(DoubletFinder) # chris-mcginnis-ucsf/DoubletFinder
library(cowplot)
library(here)
library(pheatmap)
library(dsb)
library(clustree)
library(harmony)
library(viridis)
library(pathview)
library(openxlsx)
library(gprofiler2)
library(singlecellmethods) # immunogenomics/singlecellmethods
library(batchelor)
library(knitr)
library(LaCroixColoR) # johannesbjork/LaCroixColoR
library(rmarkdown)
library(SingleCellExperiment)
library(MetBrewer) # BlakeRMills/MetBrewer
library(MoMAColors) # BlakeRMills/MoMAColors
library(SignacX)
library(djvdj)
library(markdown)
library(ComplexHeatmap)
library(RSpectra)
library(Rsamtools) # needed to download rhtslib from github not bioconductor
library(alakazam)
library(shazam)
library(tigger)
library(scoper)
library(dowser) # need GenomicAlignments and ggtree from bioconductor
library(muscat)
library(ADTnorm) # yezhengSTAT/ADTnorm
library(flowCore) 
library(flowStats)
library(EMDomics)
library(DropletUtils)
library(splitstackshape)
library(data.table)
#library(Platypus) # "alexyermanos/Platypus"
library(stringdist)
library(alluvial)
library(ggalluvial)
library(glmGamPoi)
library(scran)