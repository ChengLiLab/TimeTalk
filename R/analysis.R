#' Function to run monocle3 pipeline
#'
#' This function will run monocle3 pipeline from seurat object
#' @param tmp.seu input seurat v3 object
#' @param tmp.assay which assay to use from seurat v3 object, default is integrated
#' @param tmp.minimal_branch_len The minimal length of the diameter path for a branch to be preserved during graph pruning procedure. Default is 5.
#' @return A monocle3 cds object with trajectory information
#' @import Seurat
#' @import monocle3
#' @export
#' @examples
#' cds <- RunMonocle3(tmp.seu = tmp.seu, tmp.minimal_branch_len = 5)
RunMonocle3 <- function(tmp.seu, tmp.assay = "integrated", tmp.minimal_branch_len = 5) {
  mat <- Seurat::GetAssayData(tmp.seu, assay = "RNA", slot = "counts")
  tmp.seu.pdata <- tmp.seu[[]]
  tmp.seu.fdata <- data.frame(gene_short_name = rownames(mat))
  rownames(tmp.seu.fdata) <- rownames(mat)
  cds <- monocle3::new_cell_data_set(
    expression_data = mat,
    cell_metadata = tmp.seu.pdata,
    gene_metadata = tmp.seu.fdata
  )
  cds <- monocle3::preprocess_cds(cds = cds, num_dim = 50)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  #### load integrated umap
  tmp.cell.id <- rownames(cds@int_colData$reducedDims$UMAP)
  DefaultAssay(tmp.seu) <- tmp.assay
  tmp.umap.embed <- Seurat::Embeddings(tmp.seu, reduction = "umap")[tmp.cell.id, ]
  cds@int_colData$reducedDims$UMAP <- tmp.umap.embed
  #### clustering
  cds <- cluster_cells(cds = cds)
  ## learn trajectory
  cds <- monocle3::learn_graph(cds, close_loop = F, learn_graph_control = list(
    minimal_branch_len = tmp.minimal_branch_len,
    orthogonal_proj_tip = T
  ))
  return(cds)
}

#' Function to get variable gene in trajectory
#'
#' This function will get the variable gene along the trajectory, by Test genes for differential expression based on the low dimensional embedding and the principal graph
#' @param tmp.cds monocle3 cds object
#' @param tmp.cores the number of cores to be used while testing each gene for differential expression, default is 8
#' @param tmp.neighbor_graph String indicating what neighbor graph to use. "principal_graph" and "knn" are supported. Default is "knn", but "principal_graph" is recommended for trajectory analysis.
#' @param tmp.q_value_cutoff cutoff of q value, default is 1e-8
#' @return a named vector with Moran's I as a phenotype for the transcriptional regulatory network reconstruction, names as genes
#' @import tidyverse
#' @export
#' @examples
#' tmp.phenotype <- GetVariableGeneInTrajectory(tmp.cds = cds)
GetVariableGeneInTrajectory <- function(tmp.cds,
                                        tmp.cores = 8,
                                        tmp.neighbor_graph = "knn",
                                        tmp.q_value_cutoff = 1e-8) {
  pr_graph_test_res <- monocle3::graph_test(
    cds = tmp.cds,
    neighbor_graph = tmp.neighbor_graph,
    cores = tmp.cores
  )
  gene.use <- pr_graph_test_res %>%
    dplyr::filter(q_value < tmp.q_value_cutoff)
  pr_deg_ids <- gene.use %>%
    pull(gene_short_name) %>%
    as.vector() %>%
    as.character()
  #### use these variable gene as phenotype
  tmp.hits <- pr_graph_test_res %>%
    dplyr::filter(q_value < tmp.q_value_cutoff) %>%
    pull(gene_short_name) %>%
    as.vector() %>%
    as.character()
  tmp.phenotype <- pr_graph_test_res %>%
    dplyr::filter(q_value < tmp.q_value_cutoff) %>%
    pull(morans_test_statistic)
  names(tmp.phenotype) <- tmp.hits
  return(tmp.phenotype)
}

#' Function to get cell types specific network
#'
#' This function will genereate cell type specific network by RTN package, the TRN network reconstruction require zero-conserved imputation method ALRA implemented by Seurat
#' @param tmp.seu input Seurat v3 project
#' @param tmp.ident input celltypes
#' @param gene.TF require TF list, default is downloaded from AnimalTFDB
#' @param tmp.phenotype Moran'I, you can obtained this by run [GetVariableGeneInTrajectory]
#' @return a tna object defined by RTN package, which contains master regulator and transcription factor regulatory network information
#' @import snow
#' @import RTN
#' @import Seurat
#' @import tidyverse
#' @export
#' @examples
#' GetCellTypeSpecificTRN(tmp.seu = seu, tmp.ident = "2C_like", gene.TF = gene.TF, tmp.phenotype = tmp.phenotype)
GetCellTypeSpecificTRN <- function(tmp.seu, tmp.ident, gene.TF, tmp.phenotype) {
  cat(paste0(tmp.ident, " TRN reconsturct start:"), sep = "\n")

  ### Check CellTypes
  if (!"CellType" %in% colnames(tmp.seu[[]])) {
    stop("Please add CellType annotation in seurat object!")
  }

  object.assays <- Seurat:::FilterObjects(object = tmp.seu, classes.keep = "Assay")
  if (!"alra" %in% object.assays) {
    tmp.seu <- RunALRA(tmp.seu)
  }
  #### get the imputation data
  tmp.data.alra <- as.matrix(GetAssayData(object = tmp.seu, slot = "data", assay = "alra"))

  tmp.df <- tmp.seu[[]]
  tmp.cell.id <- tmp.df %>%
    filter(CellType == tmp.ident) %>%
    rownames()
  ### run tni pipeline
  rtni <- tni.constructor(
    expData = tmp.data.alra[, tmp.cell.id],
    regulatoryElements = gene.TF
  )
  # Please set nPermutations >= 1000
  options(cluster = snow::makeCluster(spec = 10, "SOCK"))
  rtni <- tni.permutation(rtni, nPermutations = 1000, pValueCutoff = 0.01)
  rtni <- tni.bootstrap(rtni)
  stopCluster(getOption("cluster"))
  # Compute the DPI-filtered regulatory network
  rtni <- tni.dpi.filter(rtni)
  ### run tna pipeline
  rtna <- tni2tna.preprocess(rtni, phenotype = tmp.phenotype, hits = names(tmp.phenotype))
  # Run the MRA method
  rtna <- tna.mra(rtna)

  cat(paste0(tmp.ident, " TRN reconsturct end."), sep = "\n")
  return(rtna)
}


#' Function to run TimeTalk
#'
#' This function implemented TimeTalk which utilize the pseudotime information to perform cell-cell communication in of two cell types.
#' @param tmp.cds input monocle3 cds object.
#' @param tmp.seu input Seurat v3 object.
#' @param tmp.orig.ident sample infor in seurat object.
#' @param tmp.ident.1 input sender cell type.
#' @param tmp.ident.2 input recevier cell type.
#' @param LRpairs.df LRpairs dataframe provided from CellTalkDB.
#' @param tmp.mra.res a data frame contains all the cell type specific master regulator information. To see more details, please see examples
#' @param tmp.lags 	integer specifying th order of lags to include in the auxiliary regression to perform granger causal test
#' @param tmp.winsz window size of the interpolation, default is 0.1.
#' @param tmp.cores multisession of cores, Default is 10
#' @param numPts number of desired interpolated points, default is 200.
#' @param tmp.SCC.cutoff sperman correlation cutoff to get LRpairs, default is 0.2
#' @param tmp.granger.cutoff granger test cutoff, default is 1e-2
#' @return a data frame contain all the LR pairs and related TF's
#' @import Seurat
#' @import monocle3
#' @import cellAlign
#' @import lmtest
#' @import future.apply
#' @export
#' @examples
#' # ### build all the CellType specific TRN, not run it will take too long time
#' lapply(levels(seu), function(ii) {
#'   rtna <- myGetCellTypeSpecificTRN(
#'     tmp.seu = seu,
#'     tmp.ident = ii,
#'     gene.TF = gene.TF,
#'     tmp.phenotype = tmp.phenotype
#'   )
#'   saveRDS(rtna, file = myFileName(prefix = paste0("res/R/", ii, "RTN_rtna_object"), suffix = ".rds"))
#' })

#' ### not run again!!!
#' tmp.all.de.RNA <- FindAllMarkers(object = seu,
#'                                  assay = "RNA",
#'                                  only.pos = T)
#' ####read result
#' tmp.res.list <- lapply(levels(seu),FUN = function(ii){
#'   cat(ii,sep = "\n")
#'   tmp.ident <- ii
#'   tmp.files <- list.files(path = "res/R",pattern = paste0(tmp.ident,"RTN_rtna_object"))
#'   rtna <- readRDS(file = paste0("res/R/",tmp.files))
#'   mra <- tna.get(rtna, what="mra", ntop = -1)
#'   tmp.gene.use <- tmp.all.de.RNA %>%
#'     filter(cluster == tmp.ident) %>%
#'     filter(p_val_adj < 0.05) %>%
#'     pull(gene)
#'   mra.res <- mra %>%
#'     filter(Pvalue < 0.05) %>%
#'     filter(Regulon %in% tmp.gene.use) %>%
#'     mutate(group = tmp.ident)
#'   return(mra.res)
#' })
#
#' tmp.mra.res <- Reduce(rbind,tmp.res.list)
#' saveRDS(object = tmp.mra.res,file = "res/R/B_blastoid_RTN_mra_result.rds")
#' tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
#' LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)
#' TimeTalk.result <- RunTimeTalk(tmp.cds=tmp.cds,
#'              tmp.seu=tmp.seu,
#'              tmp.ident.1="EPI",
#'              tmp.ident.2 = "PE",
#'              LRpairs.df = LRpairs.df,
#'              tmp.mra.res = tmp.mra.res,tmp.winsz = 0.1,
#'              numPts = 200,tmp.SCC.cutoff = 0.2,tmp.granger.cutoff = 1e-2)


RunTimeTalk <- function(tmp.cds,
                        tmp.seu,
                        tmp.ident.1,
                        tmp.ident.2,
                        LRpairs.df,
                        tmp.mra.res,
                        tmp.orig.ident,
                        tmp.lags = 1,
                        tmp.winsz = 0.1,
                        numPts = 200,
                        tmp.cores = 10,
                        tmp.SCC.cutoff = 0.2,
                        tmp.granger.cutoff = 1e-2) {
  cat(paste0("prepare ligand and receptor genelist"), sep = "\n")
  LRpairs <- LRpairs.df$lr_pair
  Lgenelist <- LRpairs.df$ligand_gene_symbol
  Rgenelist <- LRpairs.df$receptor_gene_symbol

  ### using data in RNA assay
  tmp.data <- GetAssayData(object = seu, slot = "data", assay = "RNA")
  gene_symbols <- rownames(tmp.data)
  l.remove <- setdiff(Lgenelist, gene_symbols)
  r.remove <- setdiff(Rgenelist, gene_symbols)
  index.remove <- c(which(Lgenelist %in% l.remove), which(Rgenelist %in% r.remove))
  LRpairs <- LRpairs[-index.remove]
  Lgenelist <- Lgenelist[-index.remove]
  Rgenelist <- Rgenelist[-index.remove]


  cat(paste0(tmp.ident.1, "-", tmp.ident.2, " start:"), sep = "\n")
  tmp.df <- data.frame(
    pseudotime = pseudotime(cds, reduction_method = "UMAP"),
    stringsAsFactors = F
  )
  tmp.df <- cbind(seu[[]], tmp.df)
  ### Check CellTypes
  if (!c("CellType") %in% colnames(tmp.df)) {
    stop("Please add CellType annotation in seurat object!")
  }

  if (!c("orig.ident") %in% colnames(tmp.df)) {
    stop("Please add CellType annotation in seurat object!")
  }

  cat(paste0("scale pseudotime"), sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.1$pseudotime
  tmp.cell.meta.1$pseudotime <- MinMaxScale(x)

  tmp.cell.meta.2 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime
  tmp.cell.meta.2$pseudotime <- MinMaxScale(x)

  tmp.mat.1 <- tmp.data[, tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[, tmp.cell.meta.2$cell_id]

  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime

  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id

  inter.tmp.mat.1 <- cellAlign::interWeights(
    expDataBatch = tmp.mat.1,
    trajCond = tmp.mat.pseudotime.1,
    winSz = tmp.winsz,
    numPts = numPts
  )
  inter.tmp.mat.2 <- cellAlign::interWeights(
    expDataBatch = tmp.mat.2,
    trajCond = tmp.mat.pseudotime.2,
    winSz = tmp.winsz,
    numPts = numPts
  )

  inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj

  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)


  plan("multisession", workers = tmp.cores)
  tmp.res.list <- future_lapply(seq_along(Lgenelist), FUN = function(ii) {
    cat(ii, sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii], ]
    y <- inter.tmp.mat.2[Rgenelist[ii], ]
    tmp.res <- cor(x, y, method = "pearson")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist, "-", Rgenelist)
  tmp.res.list.PCC <- tmp.res.list

  tmp.res.list <- future_lapply(seq_along(Lgenelist), FUN = function(ii) {
    cat(ii, sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii], ]
    y <- inter.tmp.mat.2[Rgenelist[ii], ]
    tmp.res <- cor(x, y, method = "spearman")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist, "-", Rgenelist)
  tmp.res.list.SCC <- tmp.res.list
  plan("sequential")


  tmp.res.cor <- data.frame(
    PCC = unlist(tmp.res.list),
    SCC = unlist(tmp.res.list),
    stringsAsFactors = F
  ) %>%
    rownames_to_column("LRpairs") %>%
    mutate(PCC = ifelse(is.na(PCC), 0, PCC)) %>%
    mutate(SCC = ifelse(is.na(SCC), 0, SCC)) %>%
    arrange(-SCC) %>%
    mutate(Rank = row_number())

  tmp.LR.list <- tmp.res.cor %>%
    filter(abs(SCC) > tmp.SCC.cutoff) %>%
    pull(LRpairs)

  #### read master regulator analysis
  tmp.TF.gene <- tmp.mra.res %>%
    filter(group == tmp.ident.2) %>%
    pull(Regulon)


  plan("multisession", workers = tmp.cores)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR) {
    cat(tmp.LR, sep = "\n")
    tmp.L.gene <- unlist(lapply(strsplit(tmp.LR, split = "-"), FUN = function(ii) {
      ii[1]
    }))
    tmp.R.gene <- unlist(lapply(strsplit(tmp.LR, split = "-"), FUN = function(ii) {
      ii[2]
    }))

    x <- inter.tmp.mat.1[tmp.L.gene, ]
    y <- inter.tmp.mat.2[tmp.R.gene, ]

    tmp.res.list <- lapply(tmp.TF.gene, function(tmp.TF.gene.use) {
      cat(paste0(tmp.ident.2, "_TF:", tmp.TF.gene.use), sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use, ]
      tmp.IS <- sqrt(x * y)
      tmp.res.LRtoTF.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.IS, tmp.TF.level, order = tmp.lags)
          tmp.res$`Pr(>F)`[2]
        },
        error = function(e) {
          1
        }
      )

      tmp.res.TFtoLR.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.TF.level, tmp.IS, order = tmp.lags)
          tmp.res$`Pr(>F)`[2]
        },
        error = function(e) {
          1
        }
      )
      tmp.PCC <- cor(tmp.IS, tmp.TF.level, method = "pearson")
      tmp.PCC <- ifelse(is.na(tmp.PCC), 0, tmp.PCC)
      tmp.SCC <- cor(tmp.IS, tmp.TF.level, method = "spearman")
      tmp.SCC <- ifelse(is.na(tmp.SCC), 0, tmp.SCC)
      tmp.res.df <- data.frame(
        LR = tmp.LR,
        L = tmp.L.gene,
        R = tmp.R.gene,
        TF = tmp.TF.gene.use,
        PCC = tmp.PCC,
        SCC = tmp.SCC,
        LRtoTF = tmp.res.LRtoTF.pvalue,
        TFtoLR = tmp.res.TFtoLR.pvalue,
        cell.L = tmp.ident.1,
        cell.R = tmp.ident.2,
        stringsAsFactors = F
      )
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind, tmp.res.list)
    tmp.res.df <- tmp.res.df %>%
      mutate(category = ifelse(LRtoTF < tmp.granger.cutoff | TFtoLR < tmp.granger.cutoff, "PASS", "SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  tmp.ttt.res.df <- Reduce(rbind, tmp.ttt.res)
  tmp.res.df <- tmp.ttt.res.df
  return(tmp.res.df)
}
