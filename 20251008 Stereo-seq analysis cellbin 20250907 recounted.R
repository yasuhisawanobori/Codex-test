library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(ggplot2)
library(SeuratDisk)
library(dplyr)

# Read and convert the cellbin.h5ad file----
library(reticulate)
anndata <- import("anndata", convert = FALSE)

in_h5ad <- "C:/raw data/20250906_RatThyF-1_recount_snfix/RatThyF-1/outs/analysis/Y01052GC.cellbin_1.0.adjusted.h5ad"

# Make an explicit output path in the SAME folder
out_h5ad <- sub("\\.h5ad$", ".noobsm.h5ad", in_h5ad)

# Sanity checks (shapes and gene order)
ad <- anndata$read_h5ad(in_h5ad)
ad$obsm$clear()
ad$write_h5ad(out_h5ad)

reticulate::py_to_r(ad$layers$keys())   # look for "counts" / "raw" etc.
!is.null(ad$raw)                        # some files keep raw in .raw.X

same_genes <- all(py_to_r(ad$raw$var_names) == py_to_r(ad$var_names)) # Same number/order of genes?
py_to_r(ad$X$shape); py_to_r(ad$raw$X$shape); same_genes # If same_genes is TRUE and shapes match, proceed.

# Read H5AD → SCE. If your zellkonverter supports it, skip reduced dims explicitly.
## 1) Point basilisk to a clean folder (no spaces), then verify
dir.create("C:/basilisk", showWarnings = FALSE)
Sys.setenv(BASILISK_EXTERNAL_DIR = "C:/basilisk")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("basilisk.utils", quietly = TRUE)) BiocManager::install("basilisk.utils")

# sanity checks (both should print C:/basilisk or a path inside it)
basilisk.utils::getExternalDir()
basilisk.utils::getCondaDir()  # may not exist yet, that's OK

## 2) Install Miniforge into that folder (now space-free)
basilisk.utils::installConda()

## 3) Read your file (either path works now)
library(zellkonverter)

sce <- readH5AD(out_h5ad)                  # normal

# search raw counts data
library(Matrix)
has_decimals <- function(m) {
  x <- m@x; if (length(x)==0) return(FALSE)
  any(abs(x - round(x)) > 1e-8)
}

assayNames(sce)
integral_like <- vapply(assayNames(sce),
                        function(nm) !has_decimals(assay(sce, nm)),
                        logical(1))
integral_like

# Make sure we have a 'counts' assay
print(assayNames(sce))
if (!"counts" %in% assayNames(sce)) {
  if ("X" %in% assayNames(sce)) assay(sce, "counts") <- assay(sce, "X")
}

# Convert to Seurat
CellBin_obj <- as.Seurat(sce, counts = "counts", data = NULL)
CellBin_obj  <- UpdateSeuratObject(CellBin_obj)

# Wire .raw.X into Seurat----
# ad is your AnnData object already loaded
# 1) Convert AnnData .raw.X (cells×genes) → genes×cells CSC, then to R's dgCMatrix
raw_csc <- ad$raw$X$T$tocsc()
counts  <- reticulate::py_to_r(raw_csc)   # dgCMatrix in R

# 2) Attach dimnames (genes, cells) from AnnData's .raw var/obs
genes <- reticulate::py_to_r(ad$raw$var_names$to_list())
cells <- reticulate::py_to_r(ad$obs_names$to_list())
dimnames(counts) <- list(genes, cells)


# 3) Drop into Seurat as real counts
CellBin_obj[["originalexp"]] <- CreateAssayObject(counts = counts)
DefaultAssay(CellBin_obj) <- "originalexp"

# 4) Quick sanity checks
stopifnot(!any(abs(CellBin_obj[["originalexp"]]@counts@x - 
                     round(CellBin_obj[["originalexp"]]@counts@x)) > 1e-8))
stopifnot(all(Matrix::colSums(CellBin_obj[["originalexp"]]@counts) > 0))

# Create a clean RNA assay from originalexp counts
rna_counts <- GetAssayData(CellBin_obj, assay = "originalexp", layer = "counts")
CellBin_obj[["RNA"]] <- CreateAssayObject(counts = rna_counts)
DefaultAssay(CellBin_obj) <- "RNA"

#rename features using rowData(CellBin_obj)$real_gene_name----
#Convert the assay to Assay5
assay <- CellBin_obj[["RNA"]]
assay <- as(assay, Class = "Assay5")   # needs SeuratObject v5+

#Apply your new names
feat <- as.data.frame(rowData(sce))
symvec <- as.character(feat$real_gene_name)          # <- ensure character
id2sym <- setNames(symvec, rownames(feat))

rn <- rownames(assay)

mapped <- as.character(id2sym[rn])                   # character, may contain NA
mapped[is.na(mapped)] <- ""                          # treat NA as empty

new <- ifelse(mapped != "", mapped, rn)              # use symbol when available
new <- make.unique(new)                              # avoid duplicated names

rownames(assay) <- new

#Put back and continue
CellBin_obj[["RNA"]] <- assay
DefaultAssay(CellBin_obj) <- "RNA"

#Try to locate spatial coordinates in colData----
md <- CellBin_obj@meta.data
coord_candidates <- grep("(^x$|^y$|_x$|_y$|col$|row$|imagecol|imagerow|^X$|^Y$|spatial.*x|spatial.*y)",
                         names(md), ignore.case = TRUE, value = TRUE)

#Heuristic: pick two distinct columns for x,y
pick_xy <- function(df, cands){
  # try common pairs first
  pairs <- list(c("imagecol","imagerow"),
                c("x","y"), c("X","Y"),
                c("spatial_x","spatial_y"),
                c("col","row"))
  for (p in pairs) {
    if (all(p %in% names(df))) return(df[,p])
  }
  # fallback: first two unique candidates
  df[, unique(cands)[seq_len(2)]]
}
coords <- pick_xy(md, coord_candidates)
colnames(coords) <- c("x","y")

#Register a "spatial" DimReduc so you can use DimPlot like scRNA
colnames(coords) <- paste0("SP_", seq_len(ncol(coords)))  # SP_1, SP_2
spatial_dr <- CreateDimReducObject(embeddings = as.matrix(coords),
                                   key = "SP_", assay = DefaultAssay(CellBin_obj))
CellBin_obj[["spatial"]] <- spatial_dr

# Basic scRNA-seq workflow----
CellBin_obj <- NormalizeData(CellBin_obj)
CellBin_obj <- FindVariableFeatures(CellBin_obj, nfeatures=2000)
CellBin_obj <- ScaleData(CellBin_obj, features = VariableFeatures(CellBin_obj))

# 1) Make a sketch
set.seed(1)
n_ref <- 100000  # tune 100k–200k
ref_cells <- sample(colnames(CellBin_obj), n_ref)
ref <- subset(CellBin_obj, cells = ref_cells)

# 2) Do PCA/UMAP/clustering on the sketch
ref <- ScaleData(ref, features = VariableFeatures(CellBin_obj), verbose = TRUE)
ref <- RunPCA(ref, features = VariableFeatures(CellBin_obj),
              npcs = 40, approx = TRUE, verbose = TRUE)

ref <- FindNeighbors(ref, reduction = "pca", dims = 1:40,
                     nn.method = "annoy", n.trees = 50)
ref <- FindClusters(ref, algorithm = 4, resolution = 0.3)

ref <- RunUMAP(ref, reduction = "pca", dims = 1:40,
               umap.method = "uwot", n.neighbors = 30, min.dist = 0.3,
               return.model = TRUE)

# 3A) FAST path: map the full dataset via anchors with PCA *projection*
#     (no full Scale/PCA/UMAP needed)
anchors <- FindTransferAnchors(
  reference = ref,
  query = CellBin_obj,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "pcaproject",    # <— key: project query onto ref PCs
  dims = 1:40,
  features = VariableFeatures(CellBin_obj)
)

CellBin_obj <- MapQuery(
  anchorset = anchors,
  query = CellBin_obj,
  reference = ref,
  refdata = list(cluster = Idents(ref)),  # transfers cluster labels
  reference.reduction = "pca",
  reduction.model = "umap",               # projects into ref UMAP
  verbose = TRUE
)

# Results:
# - CellBin_obj[["ref.pca"]]  : projected PCs for all cells
# - CellBin_obj[["ref.umap"]] : projected UMAP for all cells
# - CellBin_obj$predicted.cluster : transferred cluster labels

# 3B) If you prefer not to run anchors:
#     You can still project UMAP if you have PCs for all cells.
#     (e.g., if you later compute PCs, use:)
# emb_full <- uwot::umap_transform(
#   Embeddings(CellBin_obj, "pca")[, 1:40],
#   model = ref@reductions$umap@misc$model
# )
# CellBin_obj[["umap"]] <- CreateDimReducObject(embeddings = emb_full, key = "UMAP_")

#save the original (uncleaned) object
library(SeuratDisk)
SaveH5Seurat(CellBin_obj, filename = "20250906 Rat thymus stereo-seq recount/CellBin_obj.h5seurat", overwrite = TRUE)
CellBin_obj <- LoadH5Seurat("20250906 Rat thymus stereo-seq recount/CellBin_obj.h5seurat") #Error: Call assays must have either a 'counts' or 'data' slot, missing for RNA

# 2) Trim: keep counts+data (+ the PCA you just made), drop scale.data
clean_for_save <- function(obj) {
  # 0) Nuke recursion-prone slots
  obj@neighbors <- list()
  obj@graphs    <- list()
  obj@tools     <- list()
  obj@commands  <- list()
  obj@images    <- list()  # drop spatial images if any (keeps coords in reductions)
  
  # 1) Strip UMAP/other reduction models & misc
  for (rd in names(obj@reductions)) {
    if (!is.null(obj@reductions[[rd]]@misc)) {
      # specifically remove uwot model if present
      if (!is.null(obj@reductions[[rd]]@misc$model)) {
        obj@reductions[[rd]]@misc$model <- NULL
      }
      # and then clear misc entirely to be safe
      obj@reductions[[rd]]@misc <- list()
    }
  }
  
  # 2) Keep only lean layers & only the reductions you actually need
  keep_reds <- intersect(c("pca","umap","ref.pca","ref.umap","spatial"),
                         names(obj@reductions))
  obj <- DietSeurat(
    obj,
    assays    = "RNA",
    layers    = list(RNA = c("counts","data")),   # drop scale.data
    dimreducs = keep_reds,
    graphs    = character()
  )
  obj
}

CellBin_obj <- clean_for_save(CellBin_obj)

## Preferred: HDF5 (robust for big objects)
SaveH5Seurat(CellBin_obj, filename = "20250906 Rat thymus stereo-seq recount/CellBin_obj.h5seurat", overwrite = TRUE)

loadtest <- LoadH5Seurat("20250906 Rat thymus stereo-seq recount/CellBin_obj.h5seurat") # Error: Call assays must have either a 'counts' or 'data' slot, missing for RNA

#Set the leiden clusters as clusters of the Seurat object----
# 1) Make sure it's a clean factor (avoid weird types)
CellBin_obj$leiden <- as.character(CellBin_obj$leiden)
CellBin_obj$leiden <- ifelse(is.na(CellBin_obj$leiden), "NA", CellBin_obj$leiden)
CellBin_obj$leiden <- factor(CellBin_obj$leiden)

# 2) Set active identities to SAW's Leiden
Idents(CellBin_obj) <- "leiden"

CellBin_obj$leiden <- factor(CellBin_obj$leiden,
                             levels = sort(as.numeric(unique(CellBin_obj$leiden))))

CellBin_obj <- SetIdent(CellBin_obj, value = factor(CellBin_obj$leiden), levels = levels(unique(CellBin_obj$leiden)))

# 3) Quick overview
table(Idents(CellBin_obj))[1:20]   # first 20 clusters + sizes
length(levels(Idents(CellBin_obj)))  # number of Leiden clusters

DimPlot(CellBin_obj, reduction = "spatial", group.by = "leiden", label = TRUE)

# Quick overview----
FeaturePlot(CellBin_obj, reduction = "spatial", features = c("Epcam", "Sell", "Cd3e", "S1pr1","Cd4", "Cd8b"))
DotPlot(CellBin_obj, features = c("Epcam", "Sell", "Cd3e", "S1pr1","Cd4", "Cd8b"))

#Save----
options(expressions = 5e5)               # raise recursion limit as a safety net
saveRDS(CellBin_obj, "20250906 Rat thymus stereo-seq recount/CellBin_obj.rds") 
CellBin_obj <- readRDS("20250906 Rat thymus stereo-seq recount/CellBin_obj.rds")

# Markers (per cluster)----
markers <- FindAllMarkers(CellBin_obj, only.pos = TRUE, logfc.threshold = 0.25)
head(markers)

# Automatic mECA/mEFA recognition----
# pick signatures (adapt/expand as you like)
mTEC_genes  <- c("Epcam","Krt5","Krt8","Reg3a","Reg3g","Alox12e", "Fezf2")
mThymo_genes <- c("Cd69", "Tox", "Cd2", "Cd28", "Ms4a4a","Cd53","Sell","Ccr7")

# (1) module scores
CellBin_obj <- AddModuleScore(CellBin_obj, features=list(mTEC_genes),  name="mTEC_")
CellBin_obj <- AddModuleScore(CellBin_obj, features=list(mThymo_genes), name="mThymo_")

# (2) normalize scores (z within object)
md <- CellBin_obj@meta.data
md$mTEC_norm  <- as.numeric(scale(md$mTEC_1))
md$mThymo_norm <- as.numeric(scale(md$mThymo_1))

# (3) simple rule to start:
#   mEFA:   mThymo high & mTEC low
#   mECA:   mThymo high & mTEC high
#   other:  everything else
q <- function(x,p=.75) as.numeric(quantile(x,p,na.rm=TRUE))

hi_mThymo <- md$mThymo_norm >= quantile(md$mThymo_norm,.90, na.rm=TRUE)
hi_mTEC   <- md$mTEC_norm  >= quantile(md$mTEC_norm,.75, na.rm=TRUE)
lo_mTEC   <- md$mTEC_norm  <= quantile(md$mTEC_norm,.74, na.rm=TRUE)

label <- ifelse(hi_mThymo & lo_mTEC, "mEFA",
              ifelse(hi_mThymo & hi_mTEC, "mECA", "cortex"))
CellBin_obj$mRegion <- factor(label, levels=c("mEFA","mECA","cortex"))

a <-DimPlot(CellBin_obj, reduction="spatial", group.by="mRegion", pt.size=0.5,
             cells     = rownames(CellBin_obj@meta.data)[CellBin_obj$mRegion %in% c("mEFA","mECA", "cortex")]) +
       scale_y_reverse() + coord_fixed()

# (4) optional spatial smoothing (reduces salt-pepper):
# use 2D kNN average of scores
library(FNN)
k <- 200
knn <- get.knn(as.matrix(md[,c("x","y")]), k=k)$nn.index
smooth <- function(v) rowMeans(matrix(v[knn], ncol=k))
md$Thymo_smooth <- smooth(md$mThymo_norm); md$mTEC_smooth <- smooth(md$mTEC_norm)
hi_mThymo   <- md$Thymo_smooth >= q(md$Thymo_smooth,.90)
hi_mTEC   <- md$mTEC_smooth  >= q(md$mTEC_smooth,.75)
lo_mTEC   <- md$mTEC_smooth  <= q(md$mTEC_smooth,.74)
CellBin_obj$mRegion <- factor(ifelse(hi_mThymo & lo_mTEC ,"mEFA", ifelse(hi_mThymo  & hi_mTEC ,"mECA","cortex")),
                              levels=c("mEFA","mECA","cortex"))

b <- DimPlot(
  CellBin_obj,
  reduction = "spatial",
  group.by  = "mRegion",
  cells     = rownames(CellBin_obj@meta.data)[CellBin_obj$mRegion %in% c("mEFA","mECA", "cortex")],
  pt.size   = 0.5
) +
  scale_y_reverse() + coord_fixed()

c <- DimPlot(
  CellBin_obj,
  reduction = "spatial",
  group.by  = "mRegion",
  cells     = rownames(CellBin_obj@meta.data)[CellBin_obj$mRegion %in% c("mEFA","mECA", "cortex")],
  pt.size   = 2
) +
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(15000, 11500))

a|b|c

# capsule recognition----
leiden10_17 <- subset(CellBin_obj, subset = leiden %in% c("10", "17"))

# Tiny Shiny “gadget”
lasso_select <- function(df) {
  stopifnot(all(c("x","y","cell") %in% names(df)))
  library(shiny); library(plotly); library(miniUI)
  
  ui <- miniPage(
    gadgetTitleBar("Lasso select"),
    miniContentPanel(
      plotlyOutput("plt", height = "100%")
    )
  )
  
  server <- function(input, output, session) {
    output$plt <- renderPlotly({
      plot_ly(
        df, x=~x, y=~y, key=~cell, customdata=~cell,
        type="scattergl", mode="markers",
        marker=list(size=2, opacity=0.8)
      ) %>%
        layout(
          dragmode="lasso",
          xaxis=list(title="SP_1"),
          yaxis=list(title="SP_2", autorange="reversed")
        )
    })
    
    observeEvent(input$done, {
      sel <- event_data("plotly_selected")
      picked <- if (is.null(sel)) character(0) else sel$key
      stopApp(picked)
    })
    observeEvent(input$cancel, { stopApp(character(0)) })
  }
  
  runGadget(ui, server, viewer = dialogViewer("Lasso", width = 1000, height = 800))
}

coords <- Embeddings(leiden10_17[["spatial"]])[,1:2]
df <- data.frame(cell = colnames(leiden10_17), x = coords[,1], y = coords[,2])

picked_cells <- lasso_select(df)
length(picked_cells)

leiden10_17$lasso_plotly <- FALSE
leiden10_17$lasso_plotly[match(picked_cells, colnames(leiden10_17))] <- TRUE

# Visualize the handwritten result
DimPlot(leiden10_17, reduction = "spatial", group.by = "lasso_plotly")+
  scale_y_reverse() + coord_fixed()

# Save the handwritten result
saveRDS(leiden10_17, "20250813 Rat thymus stereo-seq analysis/capFb_CellBin.rds")

leiden10_17 <- readRDS("20250813 Rat thymus stereo-seq analysis/capFb_CellBin.rds")

# Map the lasso-selected capFbs back to the original object
picked_cells <- rownames(leiden10_17@meta.data[leiden10_17@meta.data$lasso_plotly == TRUE,])
CellBin_obj$capsule_lasso <- FALSE
CellBin_obj$capsule_lasso[match(picked_cells, colnames(CellBin_obj))] <- TRUE

DimPlot(CellBin_obj, reduction = "spatial", group.by = "capsule_lasso") +
  scale_y_reverse() + coord_fixed()

# Expand the selection to proximal bins
## 0) pre: CellBin_obj$capsule_lasso is TRUE for your lasso picks
coords_all <- Embeddings(CellBin_obj[["spatial"]])[, 1:2]
sel_idx    <- which(CellBin_obj$capsule_lasso)

# 1) KD-tree nearest neighbor distance from ALL → SELECTED
#    (FNN is very fast and memory-light)
library(FNN)
d_nearest <- as.vector(
  FNN::knnx.dist(data  = coords_all[sel_idx, , drop = FALSE],
                 query = coords_all,
                 k = 1)
)

## 2) choose a radius in the same units as your spatial axes
radius <- 100   # 50 µm (1 bin = 0.5 µm)

CellBin_obj$capsule_expanded <- d_nearest <= radius

## 3) visualize
DimPlot(CellBin_obj, reduction = "spatial", group.by = "capsule_expanded") +
  scale_y_reverse() + coord_fixed()

## 4) mix with mRegion
if (is.factor(CellBin_obj$mRegion)) {
  levels(CellBin_obj$mRegion) <- union(levels(CellBin_obj$mRegion), "capsular")
}

CellBin_obj$thyRegion <- as.character(CellBin_obj$mRegion)
CellBin_obj$capsule_dmin <- d_nearest

set_capsule_by_radius <- function(obj, radius, allowed_regions = NULL) {
  base <- obj$mRegion                      # fixed snapshot
  mask <- obj$capsule_dmin <= radius            # new selection
  
  # rebuild fresh labels from the snapshot, then overlay "capsular"
  new_lab <- base
  new_lab[mask] <- "capsular"
  
  # write back as factor with a stable level set
  obj$thyRegion <- factor(new_lab, levels = union("capsular", unique(base)))
  
  # keep handy metadata for audit/plotting
  obj$capsule_radius  <- radius
  obj$capsule_mask    <- mask
  
  obj
}

CellBin_obj <- set_capsule_by_radius(CellBin_obj, radius = 100)

DimPlot(CellBin_obj, reduction = "spatial", group.by = "thyRegion")+
  scale_y_reverse() + coord_fixed()

# Clustering and characterization(RCTD)----
RCTD_CellBin <-readRDS("20250906 Rat thymus stereo-seq recount/RCTD_CellBin20260114.rds")

# IDs RCTD actually used
rctd_ids <- rownames(RCTD_CellBin@spatialRNA@coords)

# ---- build a compact metadata frame from your list-style results (top1/top2)
res_list <- RCTD_CellBin@results

get_top1 <- function(x) if (length(x$cell_type_list)>=1) x$cell_type_list[1] else NA_character_
get_top2 <- function(x) if (length(x$cell_type_list)>=2) x$cell_type_list[2] else NA_character_
get_top3 <- function(x) if (length(x$cell_type_list)>=3) x$cell_type_list[3] else NA_character_
get_w    <- function(x, nm) if (!is.na(nm) && length(x$sub_weights)) unname(x$sub_weights[nm]) else NA_real_
get_conf <- function(x, nm) if (!is.na(nm) && "conf_list" %in% names(x)) unname(x$conf_list[nm]) else NA

top1 <- vapply(res_list, get_top1, character(1))
top2 <- vapply(res_list, get_top2, character(1))
top3 <- vapply(res_list, get_top3, character(1))
w1   <- mapply(function(x,nm) get_w(x,nm),   res_list, top1, USE.NAMES=FALSE)
w2   <- mapply(function(x,nm) get_w(x,nm),   res_list, top2, USE.NAMES=FALSE)
w3   <- mapply(function(x,nm) get_w(x,nm),   res_list, top3, USE.NAMES=FALSE)
c1   <- mapply(function(x,nm) get_conf(x,nm),res_list, top1, USE.NAMES=FALSE)
c2   <- mapply(function(x,nm) get_conf(x,nm),res_list, top2, USE.NAMES=FALSE)
c3   <- mapply(function(x,nm) get_conf(x,nm),res_list, top3, USE.NAMES=FALSE)
ms   <- vapply(res_list, function(x) x$min_score, numeric(1))
ca   <- vapply(res_list, function(x) x$conv_all,  logical(1))
cs   <- vapply(res_list, function(x) x$conv_sub,  logical(1))

meta_rctd <- data.frame(
  RCTD_top1      = top1,
  RCTD_top1_w    = w1,
  RCTD_top1_conf = c1,
  RCTD_top2      = top2,
  RCTD_top2_w    = w2,
  RCTD_top2_conf = c2,
  RCTD_top3      = top3,
  RCTD_top3_w    = w3,
  RCTD_top3_conf = c3,
  RCTD_min_score = ms,
  RCTD_conv_all  = ca,
  RCTD_conv_sub  = cs,
  row.names = rctd_ids,
  check.names = FALSE
)

# ---- expand to all cells (NA for ones not present in RCTD)
all_ids <- colnames(CellBin_obj)
meta_all <- meta_rctd[intersect(rownames(meta_rctd), all_ids), , drop=FALSE]   # keep only valid columns
blank <- setdiff(all_ids, rownames(meta_all))
if (length(blank)) {
  na_block <- as.data.frame(matrix(NA, nrow=length(blank), ncol=ncol(meta_rctd)),
                            row.names = blank)
  names(na_block) <- names(meta_rctd)
  meta_all <- rbind(meta_all, na_block)
}
meta_all <- meta_all[all_ids, , drop=FALSE]  # order to match object

CellBin_obj <- AddMetaData(CellBin_obj, metadata = meta_all)

# Plots
lvl <- c("Projected: DN", 
         "Projected: DP-P", "Projected: DP", "Projected: Pre DP-sig", 
         "Projected: DP-sig", "Projected: Immature CD4-1", "Projected: Immature CD4-2",
         "Projected: Mature CD4-1", "Projected: Mature CD4-2", 
         "Projected: Immature CD8-1", "Projected: Immature CD8-2", "Projected: Mature CD8-1", 
         "Projected: Mature CD8-2", "Projected: AgoSel1", "Projected: Pre-Treg",
         "Projected: Treg", "Projected: AgoSel2", 
         "capFb", "TMC1", "TMC2", "TMC3", "TMC4", "vSMC_PC", "Projected: TEC:early Pr", 
         "Projected: TEC:cTEC", "Projected: TEC:mTEC1", "Projected: TEC:mTEC-prol", 
         "Projected: TEC:mTEC2",
         "Projected: EC:capEC", "Projected: EC:aEC", 
         "Projected: EC:vEC")

CellBin_obj@meta.data$RCTD_top1 <- factor(CellBin_obj@meta.data$RCTD_top1, levels = lvl)
CellBin_obj@meta.data$RCTD_top2 <- factor(CellBin_obj@meta.data$RCTD_top2, levels = lvl)
CellBin_obj@meta.data$RCTD_top3 <- factor(CellBin_obj@meta.data$RCTD_top3, levels = lvl)

# 20 x 12 in
DimPlot(subset(CellBin_obj, subset = RCTD_top1_conf == T), reduction = "spatial",
        group.by = "RCTD_top1", split.by   = "RCTD_top1",
        pt.size    = 2, ncol = 7) + scale_y_reverse() + coord_fixed()+NoLegend()

# 20 x 12 in
DimPlot(subset(CellBin_obj, subset = RCTD_top2_conf == T), reduction = "spatial",
        group.by = "RCTD_top2", split.by   = "RCTD_top2",
        pt.size    = 1, ncol = 7) + scale_y_reverse() + coord_fixed()+NoLegend()

DimPlot(subset(CellBin_obj, subset = RCTD_top3_conf == T), reduction = "spatial",
        group.by = "RCTD_top3", split.by   = "RCTD_top3",
        pt.size    = 2, ncol = 7) + scale_y_reverse() + coord_fixed()+NoLegend()

DimPlot(subset(CellBin_obj, subset = RCTD_top1_conf == T), reduction = "spatial",
        group.by = "RCTD_top1", split.by   = "RCTD_top1",
        pt.size    = 3, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))+NoLegend()

DimPlot(subset(CellBin_obj, subset = RCTD_top2_conf == T), reduction = "spatial",
        group.by = "RCTD_top2", split.by   = "RCTD_top2",
        pt.size    = 1, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))+NoLegend()

DimPlot(subset(CellBin_obj, subset = RCTD_top3_conf == T), reduction = "spatial",
        group.by = "RCTD_top3", split.by   = "RCTD_top3",
        pt.size    = 1, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))+NoLegend()

plot_rctd_bar <- function(obj_or_md, group_col="RCTD_top1", fill_col="leiden",
                          levels=NULL, K=10, label_frac=NULL,
                          legend_title=fill_col, y_limits=NULL) {
  library(dplyr); library(ggplot2); library(tibble)
  
  # get meta
  md <- if (inherits(obj_or_md, "Seurat")) obj_or_md@meta.data else obj_or_md
  md <- as_tibble(md, rownames = "cell")
  
  # tiny helper: force atomic character (handles list-columns)
  to_chr <- functi0on(x) if (is.list(x))
    vapply(x, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
  else as.character(x)
  
  g <- to_chr(md[[group_col]])
  f <- to_chr(md[[fill_col]])
  
  # factorize x (optional ordering)
  g <- if (is.null(levels)) factor(g) else factor(g, levels = levels)
  
  # make fill numeric-ordered if it looks numeric
  if (all(grepl("^\\s*-?\\d+\\s*$", f[!is.na(f)]))) f <- factor(as.integer(f)) else f <- factor(f)
  
  md[[group_col]] <- g
  md[[fill_col]]  <- f
  
  # counts + fractions
  df <- md %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[fill_col]])) %>%
    dplyr::count(.data[[group_col]], .data[[fill_col]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  
  # labels: top-K (default) or by fraction threshold
  if (!is.null(label_frac)) {
    df <- df %>% mutate(label = if_else(frac >= label_frac, as.character(.data[[fill_col]]), NA_character_))
  } else {
    df <- df %>%
      group_by(.data[[group_col]]) %>%
      arrange(desc(n), .by_group = TRUE) %>%
      mutate(label = if_else(row_number() <= K, as.character(.data[[fill_col]]), NA_character_)) %>%
      ungroup()
  }
  
  p <- ggplot(df, aes(x=.data[[group_col]], y=n, fill=.data[[fill_col]])) +
    geom_col() +
    geom_text(aes(label=label), position=position_stack(vjust=0.5),
              size=2.8, na.rm=TRUE) +
    labs(x=group_col, y="Number of cells", fill=legend_title) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits=y_limits, expand=c(0,0))
  p
}

plot_rctd_bar(subset(CellBin_obj, subset = RCTD_top1_conf == T)@meta.data, group_col="RCTD_top1", levels=lvl, K=5, y_limits=c(0,60000))
plot_rctd_bar(subset(CellBin_obj, subset = RCTD_top2_conf == T)@meta.data, group_col="RCTD_top2", levels=lvl, K=5, y_limits=c(0,60000))
plot_rctd_bar(subset(CellBin_obj, subset = RCTD_top3_conf == T)@meta.data, group_col="RCTD_top3", levels=lvl, K=5, y_limits=c(0,60000))

plot_rctd_bar_simple <- function(
    obj_or_md,
    group_col = "RCTD_top1",
    levels    = NULL,         # desired x order
    drop_na   = TRUE,         # drop NA groups
    y_limits  = NULL,         # e.g., c(0, 200000)
    bar_fill  = "grey40",
    bar_width = 0.8,
    y_log     = FALSE         # <-- NEW: log scale option
) {
  library(dplyr); library(ggplot2); library(tibble)
  
  # 1) metadata
  md <- if (inherits(obj_or_md, "Seurat")) obj_or_md@meta.data else obj_or_md
  md <- as_tibble(md, rownames = "cell")
  
  # helper for list/factor → character
  to_chr <- function(x) {
    if (is.list(x))
      vapply(x, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
    else as.character(x)
  }
  
  g <- if (group_col %in% names(md)) to_chr(md[[group_col]]) else NA_character_
  if (drop_na) g <- g[!is.na(g)]
  
  # 2) counts
  df <- tibble(group = g) |>
    count(group, name = "n") |>
    mutate(group = if (is.null(levels)) factor(group) else factor(group, levels = levels)) |>
    arrange(group)
  
  # 3) plot
  p <- ggplot(df, aes(x = group, y = n)) +
    geom_col(fill = bar_fill, width = bar_width) +
    labs(x = group_col, y = "Number of cells") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13)
    )
  
  # 4) Y-axis scaling
  if (y_log) {
    p <- p + scale_y_log10(
      limits = y_limits,
      expand = c(0, 0),
      labels = scales::comma
    )
  } else if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, expand = c(0, 0))
  }
  
  p
}

plot_rctd_bar_simple(subset(CellBin_obj, subset = RCTD_top1_conf == T)@meta.data, group_col="RCTD_top1", levels=lvl, y_limits=c(1,300000), y_log = T)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
ggsave(filename = "pic/FigS3_2.pdf", plot = get_last_plot(), width = 7, height = 4.5)
plot_rctd_bar_simple(subset(CellBin_obj, subset = RCTD_top2_conf == T)@meta.data, group_col="RCTD_top2", levels=lvl, y_limits=c(1,300000), y_log = T)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 

plot_rctd_bar_simple(subset(CellBin_obj, subset = RCTD_top3_conf == T)@meta.data, group_col="RCTD_top3", levels=lvl, y_limits=c(1,300000), y_log = T)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 


plot_rctd_top_ratio <- function(
    obj_or_md,
    top1_col = "RCTD_top1",
    compare_col = "RCTD_top2",   # e.g. RCTD_top2 or RCTD_top3
    levels = lvl,
    drop_na = TRUE,
    ratio_label = NULL,
    y_limits = NULL,
    bar_fill = "grey40",
    bar_width = 0.8,
    y_log1_centered = T,
    filter_conf_false = T,
    top1_conf_col = paste0(top1_col, "_conf"),
    compare_conf_col = paste0(compare_col, "_conf")
) {
  library(dplyr); library(ggplot2); library(tibble)
  
  md <- if (inherits(obj_or_md, "Seurat")) obj_or_md@meta.data else obj_or_md
  md <- as_tibble(md)
  
  to_chr <- function(x) {
    if (is.list(x)) {
      vapply(x, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
    } else {
      as.character(x)
    }
  }
  
  to_lgl <- function(x, n) {
    if (is.null(x)) return(rep(NA, n))
    if (is.list(x)) {
      vapply(x, function(y) if (length(y)) as.logical(y[[1]]) else NA, logical(1))
    } else {
      as.logical(x)
    }
  }
  
  n <- nrow(md)
  top1 <- if (top1_col %in% names(md)) to_chr(md[[top1_col]]) else rep(NA_character_, n)
  topx <- if (compare_col %in% names(md)) to_chr(md[[compare_col]]) else rep(NA_character_, n)
  
  if (filter_conf_false) {
    top1_conf <- if (top1_conf_col %in% names(md)) to_lgl(md[[top1_conf_col]], n) else rep(NA, n)
    topx_conf <- if (compare_conf_col %in% names(md)) to_lgl(md[[compare_conf_col]], n) else rep(NA, n)
    
    top1[top1_conf %in% FALSE] <- NA_character_
    topx[topx_conf %in% FALSE] <- NA_character_
  }
  
  if (drop_na) {
    top1 <- top1[!is.na(top1)]
    topx <- topx[!is.na(topx)]
  }
  
  n_top1 <- tibble(group = top1) |>
    filter(!is.na(group)) |>
    count(group, name = "n_top1")
  
  n_topx <- tibble(group = topx) |>
    filter(!is.na(group)) |>
    count(group, name = "n_topx")
  
  df <- full_join(n_top1, n_topx, by = "group") |>
    mutate(
      n_top1 = coalesce(n_top1, 0L),
      n_topx = coalesce(n_topx, 0L),
      ratio = if_else(n_top1 > 0, n_topx / n_top1, NA_real_),
      ratio_plot = ratio,
      group = if (is.null(levels)) factor(group) else factor(group, levels = levels)
    ) |>
    arrange(group)
  
  y_lab <- if (is.null(ratio_label)) {
    "RCTD_topN/RCTD_top1"
  } else {
    ratio_label
  }
  
  if (y_log1_centered) {
    ratio_floor <- if (!is.null(y_limits) && length(y_limits) >= 1 && is.finite(y_limits[1]) && y_limits[1] > 0) {
      y_limits[1]
    } else {
      pos_min <- suppressWarnings(min(df$ratio[df$ratio > 0], na.rm = TRUE))
      if (is.finite(pos_min)) pos_min / 2 else 1e-4
    }
    df <- df |>
      mutate(ratio_plot = if_else(!is.na(ratio) & ratio <= 0, ratio_floor, ratio))
  }
  
  # label positions: always above bars for readability
  if (y_log1_centered) {
    df <- df |>
      mutate(
        label = if_else(is.na(ratio), NA_character_, sprintf("%.3f", ratio)),
        label_y = if_else(is.na(ratio_plot), NA_real_, ratio_plot * if_else(ratio_plot >= 1, 1.06, 1.14))
      )
  } else {
    ymax_ratio <- suppressWarnings(max(df$ratio, na.rm = TRUE))
    if (!is.finite(ymax_ratio)) ymax_ratio <- 1
    df <- df |>
      mutate(
        label = if_else(is.na(ratio), NA_character_, sprintf("%.3f", ratio)),
        label_y = if_else(is.na(ratio), NA_real_, ratio + ymax_ratio * 0.03)
      )
  }
  
  p <- ggplot(df, aes(x = group, y = if (y_log1_centered) ratio_plot else ratio)) +
    geom_col(fill = bar_fill, width = bar_width) +
    geom_label(
      aes(y = label_y, label = label),
      size = 2.8,
      label.size = 0,
      fill = "white",
      color = "black",
      label.padding = grid::unit(0.08, "lines"),
      na.rm = TRUE
    ) +
    labs(x = NULL, y = y_lab) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.margin = margin(t = 10, r = 6, b = 6, l = 6)
    ) +
    coord_cartesian(clip = "off")
  
  if (y_log1_centered) {
    p <- p + scale_y_log10(
      limits = y_limits,
      breaks = scales::breaks_log(n = 7),
      labels = scales::label_number(accuracy = 0.01),
      expand = ggplot2::expansion(mult = c(0, 0.16))
    )
  } else if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, expand = ggplot2::expansion(mult = c(0, 0.16)))
  }
  
  list(plot = p, data = df)
}

options(tibble.print_max = Inf, tibble.width = Inf)

plot_rctd_top_ratio(CellBin_obj, top1_col = "RCTD_top1",  compare_col = "RCTD_top2",
                    ratio_label = "RCTD Top2/Top1")
ggsave(filename = "pic/FigS3_3.pdf", plot = get_last_plot(), width = 7, height = 4.5)
plot_rctd_top_ratio(CellBin_obj, top1_col = "RCTD_top1",  compare_col = "RCTD_top3",
                    ratio_label = "RCTD Top3/Top1")
ggsave(filename = "pic/FigS3_4.pdf", plot = get_last_plot(), width = 7, height = 4.5)

# Extract RCTD-top1 prediction
CellBin_RCTD_obj <- subset(CellBin_obj, subset = RCTD_top1_conf == TRUE & (is.na(RCTD_top2_conf) | RCTD_top2_conf == FALSE))
CellBin_RCTD_obj@meta.data$RCTD_top1 <- factor(CellBin_RCTD_obj@meta.data$RCTD_top1, levels = lvl)

DimPlot(CellBin_RCTD_obj, reduction = "spatial",
        group.by = "RCTD_top1", split.by   = "RCTD_top1",
        pt.size    = 1, ncol = 7) + scale_y_reverse() + coord_fixed()

DimPlot(CellBin_RCTD_obj, reduction = "spatial",
        group.by = "RCTD_top1", split.by   = "RCTD_top1",
        pt.size    = 3, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))

plot_rctd_bar(CellBin_RCTD_obj@meta.data, group_col="RCTD_top1", levels=lvl, K=20)

# RCTD-predicted cell ratio in thymic regions----
plot_subset_ratio <- function(
    obj_or_md, target,
    layer_col   = "thymic_layers",
    top_cols    = c("RCTD_top1"),
    layer_lvls  = c("capsular","cortex","mECA","mEFA"),
    normalize   = c("by_layer"),        # or "by_target"
    conf        = T,
    clean_names   = TRUE,
    show_y_title  = T,               # ignored if top_title = TRUE
    top_title     = T,              # put the label at the top of the panel
    label_accuracy = 0.01,               # bar labels (percent points)
    axis_accuracy  = 0.01,               # y-axis tick precision (percent points)
    axis_breaks    = NULL,              # custom breaks (proportions), e.g. seq(0, .01, by=.002)
    axis_text_size = 12, axis_title_size = 14, label_size = 3.5, ylim_expand = 1.15
){
  normalize <- match.arg(normalize)
  
  # --- helpers
  clean <- function(x){
    x <- gsub("^\\s*Projected:\\s*","", x)
    gsub("\\s*bin\\s*ratio\\s*$","", x, ignore.case = TRUE)
  }
  to_chr <- function(x, n){
    if (is.null(x)) return(rep(NA_character_, n))
    if (is.list(x)) vapply(x, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
    else as.character(x)
  }
  
  target_vec  <- as.character(target)
  target_disp <- if (clean_names) clean(target_vec) else target_vec
  
  # --- data in
  if (inherits(obj_or_md, "Seurat")) {
    md <- obj_or_md@meta.data
  } else {
    md <- obj_or_md
  }
  md <- tibble::as_tibble(md, rownames = "bin")
  
  # harmonize columns
  n <- nrow(md)
  md[[layer_col]] <- factor(to_chr(md[[layer_col]], n), levels = layer_lvls)
  
  # per-top-column hits with optional *_conf gating
  per_col_hits <- lapply(top_cols, function(nm){
    val <- if (nm %in% names(md)) to_chr(md[[nm]], n) else rep(NA_character_, n)
    hit <- val %in% target_vec
    if (conf) {
      conf_col <- paste0(nm, "_conf")
      if (conf_col %in% names(md)) hit <- hit & (md[[conf_col]] %in% TRUE) else hit <- FALSE
    }
    hit
  })
  md$is_target <- Reduce(`|`, per_col_hits)
  
  # aggregate by layer
  by_layer <- md |>
    dplyr::filter(!is.na(.data[[layer_col]])) |>
    dplyr::group_by(.data[[layer_col]]) |>
    dplyr::summarise(
      n_total  = dplyr::n(),
      n_target = sum(is_target, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # normalization
  if (normalize == "by_layer") {
    by_layer$ratio <- with(by_layer, n_target / n_total)
    y_label <- paste0(paste(target_disp, collapse = " / "))
  } else {
    tot <- sum(by_layer$n_target, na.rm = TRUE)
    by_layer$ratio <- if (tot > 0) by_layer$n_target / tot else NA_real_
    y_label <- paste0(paste(target_disp, collapse = " / "), " distribution")
  }
  
  ymax <- ifelse(all(is.na(by_layer$ratio)), 1, max(by_layer$ratio, na.rm = TRUE) * ylim_expand)
  
  # --- plot
  p <- ggplot2::ggplot(by_layer, ggplot2::aes(x = .data[[layer_col]], y = ratio)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(ratio, accuracy = label_accuracy)),
                       vjust = -0.4, size = label_size) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = axis_accuracy),
      breaks = if (is.null(axis_breaks)) ggplot2::waiver() else axis_breaks,
      limits = c(0, ymax)
    ) +
    ggplot2::labs(x = NULL, y = if (show_y_title && !top_title) y_label else NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size)
    )
  
  if (top_title) {
    p <- p +
      ggplot2::labs(title = y_label) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(size = axis_title_size, hjust = 0.6, vjust = 1)
      )
  }
  p
}


plot_subset_ratio(CellBin_obj, c("Projected: DN"))
plot_subset_ratio(CellBin_obj, c("Projected: DP-P"))
plot_subset_ratio(CellBin_obj, c("Projected: DP")) 
plot_subset_ratio(CellBin_obj, c("Projected: Pre DP-sig"))
plot_subset_ratio(CellBin_obj, c("Projected: DP-sig"))
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD4-1"))
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD4-2"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD4-1"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD4-2"))
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD8-1"))
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD8-2"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD8-1"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD8-2"))
plot_subset_ratio(CellBin_obj, c("Projected: AgoSel1"))
plot_subset_ratio(CellBin_obj, c("Projected: Treg"))
plot_subset_ratio(CellBin_obj, c("Projected: AgoSel2"))
plot_subset_ratio(CellBin_obj, c("capFb"))
plot_subset_ratio(CellBin_obj, c("TMC1")) 
plot_subset_ratio(CellBin_obj, c("TMC2"))
plot_subset_ratio(CellBin_obj, c("TMC3"))
plot_subset_ratio(CellBin_obj, c("TMC4"))
plot_subset_ratio(CellBin_obj, c("vSMC_PC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:capEC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:vEC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:aEC"))
plot_subset_ratio(CellBin_obj, c("Projected: MEC"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:early Pr"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:cTEC"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC1"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC-prol"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC2"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mimetic"))

p1 <- plot_subset_ratio(CellBin_obj, c("Projected: DN"))+theme(axis.title.x = element_blank())  
p2 <- plot_subset_ratio(CellBin_obj, c("Projected: DP-P"))+theme(axis.title.x = element_blank())  
p3 <- plot_subset_ratio(CellBin_obj, c("Projected: DP"))+theme(axis.title.x = element_blank())  
p4 <- plot_subset_ratio(CellBin_obj, c("Projected: DP-sig"))+theme(axis.title.x = element_blank())  
p5 <- plot_subset_ratio(CellBin_obj, c("Projected: Immature CD4-1"))+theme(axis.title.x = element_blank())  
p6 <- plot_subset_ratio(CellBin_obj, c("Projected: Immature CD4-2"))+theme(axis.title.x = element_blank())  
p7 <- plot_subset_ratio(CellBin_obj, c("Projected: Mature CD4-1"))+theme(axis.title.x = element_blank()) 
p8 <- plot_subset_ratio(CellBin_obj, c("Projected: Mature CD4-2"))+theme(axis.title.x = element_blank())
p9 <- plot_subset_ratio(CellBin_obj, c("Projected: Immature CD8-1"))+theme(axis.title.x = element_blank())  
p10 <- plot_subset_ratio(CellBin_obj, c("Projected: Immature CD8-2"))+theme(axis.title.x = element_blank())  
p11 <- plot_subset_ratio(CellBin_obj, c("Projected: Mature CD8-1"))+theme(axis.title.x = element_blank())  
p12 <- plot_subset_ratio(CellBin_obj, c("Projected: Mature CD8-2"))+theme(axis.title.x = element_blank()) 
p13 <- plot_subset_ratio(CellBin_obj, c("Projected: AgoSel1"))+theme(axis.title.x = element_blank())  
p14 <- plot_subset_ratio(CellBin_obj, c("Projected: Pre-Treg"))+theme(axis.title.x = element_blank())  
p15 <- plot_subset_ratio(CellBin_obj, c("Projected: Treg"))+theme(axis.title.x = element_blank())  


(p1|p2|p3|p4)/(p5|p6|p7|p8)/(p9|p10|p11|p12)
ggsave(filename = "pic/Fig2A.pdf", plot = get_last_plot(), width = 10.5, height = 6.5)

(p13|p14|p15)
ggsave(filename = "pic/Fig2A_2.pdf", plot = get_last_plot(), width = 8, height = 2.5)

plot_subset_ratio <- function(
    obj_or_md, target,
    layer_col   = "thymic_layers",
    top_cols    = c("RCTD_top1", "RCTD_top2"),
    layer_lvls  = c("capsular","cortex","mECA","mEFA"),
    normalize   = c("by_layer"),        # or "by_target"
    conf        = T,
    clean_names   = TRUE,
    show_y_title  = T,               # ignored if top_title = TRUE
    top_title     = T,              # put the label at the top of the panel
    label_accuracy = 0.01,               # bar labels (percent points)
    axis_accuracy  = 0.01,               # y-axis tick precision (percent points)
    axis_breaks    = NULL,              # custom breaks (proportions), e.g. seq(0, .01, by=.002)
    axis_text_size = 12, axis_title_size = 14, label_size = 3.5, ylim_expand = 1.15
){
  normalize <- match.arg(normalize)
  
  # --- helpers
  clean <- function(x){
    x <- gsub("^\\s*Projected:\\s*","", x)
    gsub("\\s*bin\\s*ratio\\s*$","", x, ignore.case = TRUE)
  }
  to_chr <- function(x, n){
    if (is.null(x)) return(rep(NA_character_, n))
    if (is.list(x)) vapply(x, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
    else as.character(x)
  }
  
  target_vec  <- as.character(target)
  target_disp <- if (clean_names) clean(target_vec) else target_vec
  
  # --- data in
  if (inherits(obj_or_md, "Seurat")) {
    md <- obj_or_md@meta.data
  } else {
    md <- obj_or_md
  }
  md <- tibble::as_tibble(md, rownames = "bin")
  
  # harmonize columns
  n <- nrow(md)
  md[[layer_col]] <- factor(to_chr(md[[layer_col]], n), levels = layer_lvls)
  
  # per-top-column hits with optional *_conf gating
  per_col_hits <- lapply(top_cols, function(nm){
    val <- if (nm %in% names(md)) to_chr(md[[nm]], n) else rep(NA_character_, n)
    hit <- val %in% target_vec
    if (conf) {
      conf_col <- paste0(nm, "_conf")
      if (conf_col %in% names(md)) hit <- hit & (md[[conf_col]] %in% TRUE) else hit <- FALSE
    }
    hit
  })
  md$is_target <- Reduce(`|`, per_col_hits)
  
  # aggregate by layer
  by_layer <- md |>
    dplyr::filter(!is.na(.data[[layer_col]])) |>
    dplyr::group_by(.data[[layer_col]]) |>
    dplyr::summarise(
      n_total  = dplyr::n(),
      n_target = sum(is_target, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # normalization
  if (normalize == "by_layer") {
    by_layer$ratio <- with(by_layer, n_target / n_total)
    y_label <- paste0(paste(target_disp, collapse = " / "))
  } else {
    tot <- sum(by_layer$n_target, na.rm = TRUE)
    by_layer$ratio <- if (tot > 0) by_layer$n_target / tot else NA_real_
    y_label <- paste0(paste(target_disp, collapse = " / "), " distribution")
  }
  
  ymax <- ifelse(all(is.na(by_layer$ratio)), 1, max(by_layer$ratio, na.rm = TRUE) * ylim_expand)
  
  # --- plot
  p <- ggplot2::ggplot(by_layer, ggplot2::aes(x = .data[[layer_col]], y = ratio)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(ratio, accuracy = label_accuracy)),
                       vjust = -0.4, size = label_size) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = axis_accuracy),
      breaks = if (is.null(axis_breaks)) ggplot2::waiver() else axis_breaks,
      limits = c(0, ymax)
    ) +
    ggplot2::labs(x = NULL, y = if (show_y_title && !top_title) y_label else NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = axis_text_size, angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size)
    )
  
  if (top_title) {
    p <- p +
      ggplot2::labs(title = y_label) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(size = axis_title_size, hjust = 0.6, vjust = 1)
      )
  }
  p
}

plot_subset_ratio(CellBin_obj, c("Projected: DN"))
plot_subset_ratio(CellBin_obj, c("Projected: DP-P"))
plot_subset_ratio(CellBin_obj, c("Projected: DP")) 
plot_subset_ratio(CellBin_obj, c("Projected: DP-Sig"))
plot_subset_ratio(CellBin_obj, c("Projected: Neg. Sel."))
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD4"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD4"))
plot_subset_ratio(CellBin_obj, c("Projected: Treg")) 
plot_subset_ratio(CellBin_obj, c("Projected: Immature CD8"))
plot_subset_ratio(CellBin_obj, c("Projected: Mature CD8")) 
plot_subset_ratio(CellBin_obj, c("Projected: Mature cycling T cell"))
plot_subset_ratio(CellBin_obj, c("capFb"))
plot_subset_ratio(CellBin_obj, c("TMC1")) 
plot_subset_ratio(CellBin_obj, c("TMC2"))
plot_subset_ratio(CellBin_obj, c("TMC3"))
plot_subset_ratio(CellBin_obj, c("TMC4"))
plot_subset_ratio(CellBin_obj, c("vSMC_PC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:capEC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:vEC"))
plot_subset_ratio(CellBin_obj, c("Projected: EC:aEC"))
plot_subset_ratio(CellBin_obj, c("Projected: MEC"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:early Pr"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:cTEC"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC1"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC-prol"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mTEC2"))
plot_subset_ratio(CellBin_obj, c("Projected: TEC:mimetic"))

p1 <- plot_subset_ratio(CellBin_obj, c("capFb"))+theme(axis.title.x = element_blank())
p2 <- plot_subset_ratio(CellBin_obj, c("TMC1"))+theme(axis.title.x = element_blank())    
p3 <- plot_subset_ratio(CellBin_obj, c("TMC2"))+theme(axis.title.x = element_blank())
p4 <- plot_subset_ratio(CellBin_obj, c("TMC3"))+theme(axis.title.x = element_blank())
p5 <- plot_subset_ratio(CellBin_obj, c("TMC4"))+theme(axis.title.x = element_blank())      
p6 <- plot_subset_ratio(CellBin_obj, c("vSMC_PC"))+theme(axis.title.x = element_blank())      

(p1|p2|p3)/(p4|p5|p6) 

# Thymocyte-GSEA (Based on AverageExpression)----
agg_means <- AverageExpression(CellBin_obj, group.by="thyRegion", assays="RNA",  slot="data")

ranks <- (agg_means$RNA[,"mEFA"])-(agg_means$RNA[,"mECA"])
names(ranks) <- rownames(agg_means$RNA)
ranks <- sort(ranks, decreasing=TRUE)

# make reference marker table for subtraction
RatThymicRef <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicRef.rds")
RatThymicRef <- UpdateSeuratObject(RatThymicRef)
SimpleRatThymicRef <- RenameIdents(RatThymicRef, "DP1"="DP", "DP2"="DP", "Early SP"="SP", "Mature SP"="SP",
                                   "DP-P1"= "DP-P", "DP-P2"= "DP-P", "DP-P3"= "DP-P",
                                   "TMC1"="TMC1-2", "TMC2"="TMC1-2", "TMC3"="TMC3-4","TMC4"="TMC3-4",
                                   "Endo-1"="Endo", "Endo-2"="Endo",
                                   "TEC-a"="cTEC", "TEC-b"="cTEC", "TEC-c"="cTEC", "TEC-d"="mTEC", "TEC-e"="mimetic cells")

DimPlot(SimpleRatThymicRef, reduction = "umap.PCA", label = T, label.box = T, repel = T) +NoLegend()

RatThymicMarker <- FindAllMarkers(JoinLayers(SimpleRatThymicRef), logfc.threshold = 1, min.pct = 0.1, only.pos = T)
RatThymicMarker <- RatThymicMarker %>% filter(avg_log2FC > 1, pct.1 > 0.1) %>% arrange(cluster, desc(avg_log2FC))

ranks_TMC_mTECdeleted <- ranks[!(names(ranks) %in% unique(na.omit(RatThymicMarker$gene[RatThymicMarker$cluster %in% c("mTEC", "TMC3-4")])))]

# mECA thymocyte vs mEFA thymocyte scores in reference scRNA-seq
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

RatThymocyteSimple.obj <- subset(RatThymocyte.obj, subset = projected_label %in% c("Projected: Immature CD4-1", "Projected: Immature CD4-2", 
                                                                                   "Projected: Mature CD4-1", "Projected: Mature CD4-2",  
                                                                                   "Projected: Immature CD8-1", "Projected: Immature CD8-2",
                                                                                   "Projected: Mature CD8-1", "Projected: Mature CD8-2",
                                                                                   "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg"
))
RatThymocyteSimple.obj <- subset(RatThymocyteSimple.obj, subset = !is.na(projected_label))

#library(plyr)

#RatThymocyteSimple.obj$projected_label <- plyr::revalue(RatThymocyteSimple.obj$projected_label,
#                                                        c("Projected: intCD4"="Projected: Immature CD4-2",
#                                                          "Projected: Mature cycling T"="Projected: Mature CD8-2",
#                                                          "Projected: Pre-Treg"="Projected: Treg", "Projected: NKT"="Projected: AgoSel2"
#                                                          ))

mEFAthymoMarker <- names(ranks_TMC_mTECdeleted)[order(ranks_TMC_mTECdeleted, decreasing = TRUE)[seq_len(min(100, length(ranks_TMC_mTECdeleted)))]]
mECAthymoMarker <- names(ranks_TMC_mTECdeleted)[order(ranks_TMC_mTECdeleted, decreasing = F)[seq_len(min(100, length(ranks_TMC_mTECdeleted)))]]

RatThymocyteSimple.obj <- AddModuleScore(
  object   = RatThymocyteSimple.obj,
  features = list(mEFAthymoMarker),
  name     = "mEFA thymocyte score")

RatThymocyteSimple.obj <- AddModuleScore(
  object   = RatThymocyteSimple.obj,
  features = list(mECAthymoMarker),
  name     = "mECA thymocyte score")

colnames(RatThymocyteSimple.obj@meta.data)[colnames(RatThymocyteSimple.obj@meta.data)=="mECA thymocyte score1"] <- "mECA thymocyte score"
colnames(RatThymocyteSimple.obj@meta.data)[colnames(RatThymocyteSimple.obj@meta.data)=="mEFA thymocyte score1"] <- "mEFA thymocyte score"

FeaturePlot(RatThymocyteSimple.obj, features = c("mECA thymocyte score", "mEFA thymocyte score"))

# Fig. 2B, 5.2 x 5 in
RatThymocyteSimple.obj$projected_label <- factor(RatThymocyteSimple.obj$projected_label, levels = c("Projected: Immature CD4-1", "Projected: Immature CD4-2", 
                                                                                                    "Projected: Mature CD4-1", "Projected: Mature CD4-2",  
                                                                                                    "Projected: Immature CD8-1", "Projected: Immature CD8-2",
                                                                                                    "Projected: Mature CD8-1", "Projected: Mature CD8-2",
                                                                                                    "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg"
                                                                                                         ))

DotPlot(RatThymocyteSimple.obj, features = c("mECA thymocyte score", "mEFA thymocyte score"), group.by = "projected_label")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+coord_flip()
ggsave(filename = "pic/Fig2C.pdf", plot = get_last_plot(), width = 6.9, height = 2.7)

# GSEA
library(clusterProfiler);library(org.Rn.eg.db); library(enrichplot)
EFAvsECAthymoGSEA <- gseGO(geneList = ranks_TMC_mTECdeleted, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL"
                           #, pAdjustMethod="none"
                           )

gseaplot(EFAvsECAthymoGSEA, "GO:0005840", color.line="Black")
ggsave(filename = "pic/Fig2E_2.pdf", plot = get_last_plot(), width = 5.9, height = 4.2)

EFAthymoGSEA <- EFAvsECAthymoGSEA
EFAthymoGSEA@result <- EFAthymoGSEA @result[EFAthymoGSEA @result$NES > 0,]
ECAthymoGSEA <- EFAvsECAthymoGSEA
ECAthymoGSEA@result <- ECAthymoGSEA@result[ECAthymoGSEA@result$NES < 0,]

View(EFAthymoGSEA@result)

set.seed(123)
emapplot(pairwise_termsim(EFAthymoGSEA), showCategory = 1200, group=T, nCluster = 5)

treeplot(pairwise_termsim(EFAthymoGSEA), showCategory = 60, color = "p.adjust", label_format = 40,
         title = "GO GSEA (↑ in mEFA)") 

treeplot(pairwise_termsim(ECAthymoGSEA), showCategory = 40, color = "p.adjust", label_format = 40,
         title = "GO GSEA (↑ in mECA)")

dotplot(EFAthymoGSEA, showCategory = 20)+ggtitle("mEFA enriched processes")

# GSEA between rat mature and immature thymocytes
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")
RatThymocyteIN.obj <- Seurat::SetIdent(RatThymocyte.obj, value=RatThymocyte.obj$projected_label)
RatThymocyteIN.obj <- RenameIdents(RatThymocyteIN.obj, 
                                   "Projected: Mature CD4-1"="Mat", "Projected: Mature CD4-2"="Mat",
                                   "Projected: Mature CD8-1"="Mat", 
                                   #"Projected: Pre-Treg"="Mat", 
                                   "Projected: Treg"="Mat",
                                   #"Projected: DP-sig"="Imm", 
                                   "Projected: Immature CD4-1"="Imm",
                                   "Projected: Immature CD8-1"="Imm", "Projected: AgoSel1"="Imm")
RatThymoMat <- AverageExpression(RatThymocyteIN.obj, assays="RNA",  slot="data")

RatThy_MatureImmature_markers <- (RatThymoMat$RNA[,"Mat"])-(RatThymoMat$RNA[,"Imm"])
names(RatThy_MatureImmature_markers) <- rownames(RatThymoMat$RNA)
RatThy_MatureImmature_markers <- sort(RatThy_MatureImmature_markers, decreasing=TRUE)

RatThy_MatureImmature_GESA <- gseGO(geneList =RatThy_MatureImmature_markers, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL"
                                    #, pAdjustMethod="none"
                                    )


RatThy_Mature_GESA <- RatThy_MatureImmature_GESA; RatThy_Immature_GESA <- RatThy_MatureImmature_GESA
RatThy_Mature_GESA@result <- RatThy_Mature_GESA@result[RatThy_Mature_GESA@result$enrichmentScore > 0,]
RatThy_Immature_GESA@result <- RatThy_Immature_GESA@result[RatThy_Immature_GESA@result$enrichmentScore < 0,]

emapplot(pairwise_termsim(RatThy_Mature_GESA), showCategory=120)
View(RatThy_Mature_GESA@result)
treeplot(pairwise_termsim(RatThy_Mature_GESA)) | treeplot(pairwise_termsim(RatThy_Immature_GESA))
gseaplot2(RatThy_MatureImmature_GESA, geneSetID = c("GO:0005840", "GO:0022626", "GO:0003735"))

# Comparison of mEFA/ECA GSEA and Mature/Immature CD4 GSEA
commonGSEA <- EFAthymoGSEA
common_ids <- intersect(EFAthymoGSEA@result$ID, RatThy_Mature_GESA@result$ID)
commonGSEA@result <- 
  EFAthymoGSEA@result %>% filter(ID %in% common_ids)

#Fig. S4, 10x6 in
set.seed(123)
emapplot(pairwise_termsim(commonGSEA), showCategory = 50, group=T)
ggsave(filename = "pic/Fig2E_1.pdf", plot = get_last_plot(), width = 7.7, height = 4.8)

common2GSEA <- RatThy_Mature_GESA
common_ids <- intersect(RatThy_Mature_GESA@result$ID, EFAthymoGSEA@result$ID)
common2GSEA@result <- 
  RatThy_Mature_GESA@result %>% filter(ID %in% common_ids)
emapplot(pairwise_termsim(common2GSEA), showCategory = 50, group=T, nCluster = 5)

EFAthymoUniqueGSEA <- EFAthymoGSEA
EFAonly_ids <- setdiff(EFAthymoGSEA@result$ID, RatThy_Mature_GESA@result$ID)
EFAthymoUniqueGSEA@result <- EFAthymoGSEA@result %>% filter(ID %in% EFAonly_ids)
emapplot(pairwise_termsim(EFAthymoUniqueGSEA), showCategory = 50, group=T, nCluster = 5)

Thy_MatUnique_GESA <- RatThy_Mature_GESA
ThyMatonly_ids <- setdiff(RatThy_Mature_GESA@result$ID, EFAthymoGSEA@result$ID)
Thy_MatUnique_GESA@result <- RatThy_Mature_GESA@result %>% filter(ID %in% ThyMatonly_ids)
emapplot(pairwise_termsim(Thy_MatUnique_GESA), showCategory = 50, group=T, nCluster = 5)

compare_gsea_overlap <- function(gsea1, gsea2,
                                 label1 = "Set1", label2 = "Set2",
                                 bin_width = 2,
                                 id_priority = c("ID","Description"),
                                 alternative = "less",   # test if common ranks are smaller
                                 top_k = NULL            # optional hypergeometric on top-k
) {
  stopifnot(requireNamespace("dplyr"),
            requireNamespace("ggplot2"),
            requireNamespace("tidyr"))
  
  library(dplyr)
  
  # 1) Rank by p.adjust (lowest first); tie-break by NES (higher first)
  df1 <- as.data.frame(gsea1) %>%
    filter(!is.na(p.adjust)) %>%
    arrange(p.adjust, dplyr::desc(dplyr::coalesce(NES, 0))) %>%
    mutate(p.adjust_rank = dplyr::row_number(),
           .label = label1)
  
  df2 <- as.data.frame(gsea2) %>%
    filter(!is.na(p.adjust)) %>%
    arrange(p.adjust, dplyr::desc(dplyr::coalesce(NES, 0))) %>%
    mutate(p.adjust_rank = dplyr::row_number(),
           .label = label2)
  
  # 2) Choose join key: try ID, then Description
  by_field <- NULL
  for (k in id_priority) {
    if (k %in% names(df1) && k %in% names(df2)) {
      common_keys <- intersect(df1[[k]], df2[[k]])
      if (length(common_keys) > 0) { by_field <- k; break }
    }
  }
  if (is.null(by_field)) stop("Neither ID nor Description overlaps between the two results.")
  
  efa_df      <- df1
  efa_common  <- df1 %>% filter(.data[[by_field]] %in% common_keys)
  
  # 3) Build left-closed bins: 1–5, 6–10, ...
  N <- nrow(efa_df)
  starts <- seq(1, N, by = bin_width)
  ends   <- pmin(starts + bin_width - 1, N)
  breaks <- c(starts, max(ends) + 1)
  labels <- paste0(starts, "–", ends)
  
  efa_df <- efa_df %>%
    mutate(rank_bin = cut(p.adjust_rank, breaks = breaks, right = FALSE,
                          labels = labels, ordered_result = TRUE))
  
  efa_common <- efa_common %>%
    mutate(rank_bin = cut(p.adjust_rank, breaks = breaks, right = FALSE,
                          labels = labels, ordered_result = TRUE))
  
  # 4) Counts & plot data
  all_counts <- efa_df %>% count(rank_bin, name = "n_all")
  common_counts <- efa_common %>% count(rank_bin, name = "n_common")
  
  plot_df <- all_counts %>%
    full_join(common_counts, by = "rank_bin") %>%
    mutate(across(starts_with("n_"), ~tidyr::replace_na(.x, 0L)))
  
  # 5) Histogram (all = light, common = solid)
  p_hist <- ggplot2::ggplot(plot_df, ggplot2::aes(x = rank_bin)) +
    ggplot2::geom_col(ggplot2::aes(y = n_all), width = 0.9, alpha = 0.25) +
    ggplot2::geom_col(ggplot2::aes(y = n_common), width = 0.9) +
    ggplot2::labs(x = paste(label1, "p.adjust rank"),
                  y = "Number of pathways",
                  title = sprintf("Common pathways (%s ∩ %s) vs. all %s",
                                  label1, label2, label1)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # 6) Wilcoxon rank-sum: common vs other on raw ranks
  efa_df <- efa_df %>%
    mutate(is_common = .data[[by_field]] %in% common_keys)
  
  ranks_common <- efa_df$p.adjust_rank[efa_df$is_common]
  ranks_other  <- efa_df$p.adjust_rank[!efa_df$is_common]
  
  if (length(ranks_common) == 0 || length(ranks_other) == 0) {
    wilc <- list(p.value = NA_real_, statistic = NA, method = "wilcox.test", alternative = alternative)
    auc  <- NA_real_
  } else {
    wilc <- stats::wilcox.test(ranks_common, ranks_other, alternative = alternative, exact = FALSE)
    U <- as.numeric(wilc$statistic)
    m <- length(ranks_common); n <- length(ranks_other)
    auc <- U / (m * n)  # P(common < other)
  }
  
  # Optional: top-k over-representation of common pathways (hypergeometric)
  hyper <- NULL
  if (!is.null(top_k)) {
    k <- min(top_k, nrow(efa_df))
    top_set <- efa_df %>% slice_head(n = k)
    K <- sum(efa_df$is_common)
    Ntot <- nrow(efa_df)
    x <- sum(top_set$is_common)
    p_hyper <- phyper(q = x - 1, m = K, n = Ntot - K, k = k, lower.tail = FALSE)
    hyper <- list(top_k = k, overlap_in_top_k = x, total_common = K,
                  universe = Ntot, p_value = p_hyper)
  }
  
  # Return everything useful
  list(
    by_field   = by_field,
    common_keys = common_keys,
    plot_data  = plot_df,
    plot       = p_hist,
    wilcoxon   = list(p_value = wilc$p.value,
                      statistic = as.numeric(wilc$statistic),
                      alternative = alternative,
                      AUC = auc,
                      n_common = length(ranks_common),
                      n_other = length(ranks_other)),
    hypergeom  = hyper,
    ranked_tables = list(
      set1 = efa_df %>% select(all_of(by_field), p.adjust, NES, p.adjust_rank, rank_bin, is_common),
      set2 = df2   %>% select(all_of(by_field), p.adjust, NES, p.adjust_rank)
    )
  )
}

out <- compare_gsea_overlap(EFAthymoGSEA, RatThy_Mature_GESA, label1 = "mEFA", label2 = "Mature T cell", top_k = 20,
                            id_priority = "Description")
out$plot
out$wilcoxon

out2 <- compare_gsea_overlap(RatThy_Mature_GESA, EFAthymoGSEA,  label1 = "Mature T cells", label2 = "EFA", top_k = 20,
                            id_priority = "Description")
out2$plot
out2$wilcoxon


# Maturation/apoptosis gene expression comparison between mECA vs mEFA----
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

RatThymocyteSP.obj <- subset(RatThymocyte.obj, subset=projected_label %in%
                               c("Projected: DP-sig",
                                 "Projected: Immature CD4-1",  "Projected: Immature CD4-2",
                                 "Projected: Mature CD4-1", "Projected: Mature CD4-2",
                                 "Projected: AgoSel1","Projected: Treg", "Projected: AgoSel2",
                                 "Projected: Immature CD8-1","Projected: Immature CD8-2",
                                 "Projected: Mature CD8-1", "Projected: Mature CD8-2"))

ImmatureCD4Markers <- FindMarkers(SetIdent(RatThymocyteSP.obj, value = RatThymocyteSP.obj$projected_label),
                                ident.1="Projected: Immature CD4-1", ident.2 = "Projected: Mature CD4-2",
                                min.pct = 0.2, logfc.threshold = 1)
ImmatureCD4Markers <- ImmatureCD4Markers[order(-ImmatureCD4Markers$avg_log2FC), ]
c(rownames(head(ImmatureCD4Markers)), rownames(tail(ImmatureCD4Markers))) 


ImmatureCD8Markers <- FindMarkers(SetIdent(RatThymocyteSP.obj, value = RatThymocyteSP.obj$projected_label),
                                ident.1="Projected: Immature CD8-1", ident.2 = "Projected: Mature CD8-1",
                                min.pct = 0.2, logfc.threshold = 1)
ImmatureCD8Markers <- ImmatureCD8Markers[order(-ImmatureCD8Markers$avg_log2FC), ]
c(rownames(head(ImmatureCD8Markers)), rownames(tail(ImmatureCD8Markers)))

AgoSel1Markers <- FindMarkers(SetIdent(RatThymocyteSP.obj, value = RatThymocyteSP.obj$projected_label),
                              ident.1 = "Projected: AgoSel1", ident.2 = "Projected: Treg", min.pct = 0.2, logfc.threshold = 1)
AgoSel1Markers <- AgoSel1Markers[order(-AgoSel1Markers$avg_log2FC), ]
c(rownames(head(AgoSel1Markers)), rownames(tail(AgoSel1Markers)))

DotPlot(subset(RatThymocyte.obj, subset=projected_label %in%
               c(
                 "Projected: Immature CD4-1", "Projected: Mature CD4-2",
                 "Projected: Immature CD8-1", "Projected: Mature CD8-1",
               "Projected: AgoSel1", "Projected: Treg")),
        group.by = "projected_label", 
        features = unique(c(rownames(head(ImmatureCD4Markers)), rownames(tail(ImmatureCD4Markers)),
                          rownames(head(ImmatureCD8Markers)), rownames(tail(ImmatureCD8Markers)),
                          rownames(head(AgoSel1Markers)), rownames(tail(AgoSel1Markers))), fromLast=T)
                     ) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

CellBin_exp <- AverageExpression(subset(CellBin_obj, subset=thyRegion %in% c("mEFA", "mECA")), group.by = "thyRegion")

library(pheatmap)
pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(head(ImmatureCD4Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_1.pdf", width=6, height=1.6)

pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(head(ImmatureCD8Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_2.pdf", width=6, height=1.6)

pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(head(AgoSel1Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_3.pdf", width=6, height=1.6)

pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(tail(ImmatureCD4Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_4.pdf", width=6, height=1.6)

pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(tail(ImmatureCD8Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_5.pdf", width=6, height=1.6)

pheatmap(
  t(as.matrix(CellBin_exp$RNA[c(rownames(tail(AgoSel1Markers, 30))
  ), c("mEFA", "mECA"), drop = FALSE])),
  scale = "column", cluster_rows = F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, cellwidth = 13,
  cellheight = 13, legend = F,labels_row = "")
grDevices::dev.copy2pdf(file="pic/Fig2D_6.pdf", width=6, height=1.6)

DotPlot(subset(subset(CellBin_obj, subset=RCTD_top1 %in%
                      c("Projected: DP-sig",
                        "Projected: Immature CD4-1",  "Projected: Immature CD4-2",
                        "Projected: Mature CD4-1", "Projected: Mature CD4-2",
                        "Projected: AgoSel1","Projected: Treg", "Projected: AgoSel2",
                        "Projected: Immature CD8-1","Projected: Immature CD8-2",
                        "Projected: Mature CD8-1", "Projected: Mature CD8-2")),
               subset=thyRegion %in% c("mECA", "mEFA")),
        group.by = "RCTD_top1",
        split.by = "thyRegion", cols = c("blue", "blue") ,
        features = c("Tox2", "Hivep3", "Itm2a", "Nab2", "Srgn", "Tox", "Nr4a1", "Cd28", "Rgs2", "Cd69", "Raph1", "Cd2", #early SP
                     "Klrk1", "Ms4a4c", "Calhm6", "Cd79a", "Sell", "Ccnd2", "Cd40lg", "Rasa3", "Ms4a4a", "Slc4a4", "S1pr1", "Cd96", "Ccr7", #Mature SP
                     "Foxp3", "Il2ra", #Treg
                     "Pik3ca", "Ras"
        ), scale.max = 7) + RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

early_apoptosis <- c(
  # TCRhi → deletion signaling
  "Nr4a1","Nr4a3","Egr2","Egr3","Dusp1","Dusp5","Dusp6",
  # BH3-only initiators
  "Bcl2l11","Bbc3","Pmaip1","Bad","Bmf","Hrk",
  # Death-receptor branch (auxiliary in thymocytes)
  "Fas","Faslg","Fadd","Tradd","Casp8","Tnfrsf1a","Ripk1",
  # Pro-survival (often down in cells fated for deletion)
  "Bcl2","Bcl2l1","Mcl1"
)

# Later phase (commitment/execution: post-MOMP → caspases/DNA fragmentation)
later_apoptosis <- c(
  # MOMP & apoptosome
  "Bax","Bak1","Cycs","Apaf1","Casp9",
  # Executioner caspases
  "Casp3","Casp7",
  # Mitochondrial pro-apoptotic cofactors
  "Diablo","Htra2","Aifm1","Endog"
)


early_apoptosis_r <- convert_orthologs(early_apoptosis,input_species = "mouse", output_species = "rat",
                                       drop_nonorths = FALSE, gene_input = col, method = "gprofiler", gene_output = "columns",
                                       non121_strategy = 5)$ortholog_gene
early_apoptosis_r <- early_apoptosis_r[!early_apoptosis_r %in% "N/A"]

later_apoptosis_r <- convert_orthologs(later_apoptosis,input_species = "mouse", output_species = "rat",
                          drop_nonorths = FALSE, gene_input = col, method = "gprofiler", gene_output = "columns",
                          non121_strategy = 5)$ortholog_gene
later_apoptosis_r <- later_apoptosis_r[!later_apoptosis_r %in% "N/A"]

DotPlot(subset(RatThymicRef, subset=consolidated_cluster_s
               ==c("Projected: DP-Sig", "Projected: Neg. Sel.", 
                   "Projected: Immature CD4", "Projected: Mature CD4", "Projected: Treg", 
                   "Projected: Immature CD8", "Projected: Mature CD8", "Projected: Mature cycling T cell")),
        group.by = "consolidated_cluster_s", 
        features = c(early_apoptosis_r, later_apoptosis_r )) + 
  RotatedAxis() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
