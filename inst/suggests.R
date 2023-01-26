
# Bioconductor
bioc <- c('BiocVersion', 'Biobase', 'GenomicRanges', 'IRanges', 'TwoSampleMR',
          'VariantAnnotation', 'biomaRt', 'coloc', 'gwasvcf', 'ieugwasr',
          'openxlsx', 'phenoscanner', 'regioneR', 'reticulate', 'rtracklayer',
          'seqminer', 'tourr', 'arrayQualityMetrics', 'rGREAT', 'karyoploteR',
          'ComplexHeatmap', 'DESeq2', 'EnsDb.Hsapiens.v86',
          'FlowSorted.DLPFC.450k', 'GeneNet', 'RMariaDB', 'Rgraphviz',
          'TxDb.Hsapiens.UCSC.hg38.knownGene', 'bladderbatch', 'clusterProfiler',
          'ensembldb', 'fdrtool', 'graph', 'graphite', 'heatmaply', 'minfi',
          'org.Hs.eg.db', 'quantro', 'recount3', 'sva', 'snpStats')

bioc_dep <- c(‘formatR’, ‘gridGraphics’, ‘tweenr’, ‘polyclip’, ‘purrr’, ‘lambda.r’, ‘futile.options’, ‘fastmatch’, ‘ggfun’,
            ‘ggplotify’, ‘patchwork’, ‘ggforce’, ‘ggrepel’, ‘tidygraph’, ‘graphlayouts’, ‘tidytree’, ‘treeio’, ‘rngtools’,
            ‘GenomeInfoDbData’, ‘bitops’, ‘KEGGREST’, ‘dbplyr’, ‘filelock’, ‘affyio’, ‘gcrma’, ‘hexbin’, ‘BeadDataPackR’,
            ‘annotate’, ‘systemfonts’, ‘dichromat’, ‘futile.logger’, ‘snow’, ‘BH’, ‘HDO.db’, ‘fgsea’, ‘aplot’, ‘ggnewscale’,
            ‘ggraph’, ‘scatterpie’, ‘shadowtext’, ‘ggtree’, ‘doRNG’, ‘multtest’, ‘scrime’, ‘base64’, ‘sparseMatrixStats’,
            ‘rhdf5’, ‘rhdf5filters’, ‘Rhdf5lib’, ‘R.oo’, ‘R.methodsS3’, ‘BiocGenerics’, ‘S4Vectors’, ‘GenomeInfoDb’,
            ‘XVector’, ‘MatrixGenerics’, ‘SummarizedExperiment’, ‘Rsamtools’, ‘zlibbioc’, ‘Biostrings’, ‘AnnotationDbi’,
            ‘BSgenome’, ‘GenomicFeatures’, ‘Rhtslib’, ‘XML’, ‘BiocFileCache’, ‘RCurl’, ‘GenomicAlignments’, ‘BiocIO’,
            ‘restfulr’, ‘affy’, ‘affyPLM’, ‘beadarray’, ‘genefilter’, ‘gridSVG’, ‘hwriter’, ‘limma’, ‘setRNG’, ‘vsn’,
            ‘svglite’, ‘rjson’, ‘GetoptLong’, ‘DT’, ‘GO.db’, ‘TxDb.Hsapiens.UCSC.hg19.knownGene’, ‘doParallel’, ‘biovizBase’,
            ‘bezier’, ‘bamsignals’, ‘clue’, ‘BiocParallel’, ‘locfit’, ‘DOSE’, ‘enrichplot’, ‘GOSemSim’, ‘gson’, ‘qvalue’,
            ‘yulab.utils’, ‘AnnotationFilter’, ‘RSQLite’, ‘ProtGenerics’, ‘bumphunter’, ‘beanplot’, ‘nor1mix’, ‘siggenes’,
            ‘preprocessCore’, ‘illuminaio’, ‘DelayedMatrixStats’, ‘mclust’, ‘GEOquery’, ‘DelayedArray’, ‘HDF5Array’,
            ‘R.utils’, ‘sessioninfo’, ‘edgeR’)

install.packages("BiocManager")
BiocManager::install(c(bioc,bioc_dep))

#CRAN
cran_dep <- c(‘iterators’, ‘irlba’, ‘foreach’, ‘mixsqp’, ‘brew’, ‘SPAtest’, ‘RSpectra’, ‘gsl’, ‘bitops’, ‘rex’, ‘timechange’,
              ‘palmerpenguins’, ‘TSP’, ‘qap’, ‘gclus’, ‘ca’, ‘registry’, ‘caTools’, ‘susieR’, ‘zip’, ‘roxygen2’, ‘RcppTOML’,
              ‘here’, ‘SKAT’, ‘TeachingDemos’, ‘ash’, ‘energy’, ‘gifski’, ‘geozoo’, ‘covr’, ‘longitudinal’, ‘blob’, ‘DBI’,
              ‘lubridate’, ‘plogr’, ‘DBItest’, ‘dendextend’, ‘reshape2’, ‘seriation’, ‘webshot’, ‘assertthat’, ‘egg’, ‘gplots’,
              'gmp', ‘glmnet’)

install.packages(c("remotes",cran_dep), depend=TRUE)

# GitHub
remotes::install_github("jinghuazhao/pQTLdata")
remotes::install_github("mrcieu/ieugwasr")
remotes::install_github("mrcieu/gwasvcf")
remotes::install_github("mrcieu/TwoSampleMR")
remotes::install_github("phenoscanner/phenoscanner")
