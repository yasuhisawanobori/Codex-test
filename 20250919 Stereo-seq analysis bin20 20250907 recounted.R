library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(ggplot2)
library(SeuratDisk)

# Read and convert the bin20.h5ad file----
library(reticulate)
anndata <- import("anndata", convert = FALSE)

in_h5ad <- "C:/raw data/20250906_RatThyF-1_recount_snfix/RatThyF-1/outs/analysis/Y01052GC.bin20_1.0.h5ad"

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
bin20_obj <- as.Seurat(sce, counts = "counts", data = NULL)
bin20_obj  <- UpdateSeuratObject(bin20_obj)

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
bin20_obj[["originalexp"]] <- CreateAssayObject(counts = counts)
DefaultAssay(bin20_obj) <- "originalexp"

# 4) Quick sanity checks
stopifnot(!any(abs(bin20_obj[["originalexp"]]@counts@x - 
                     round(bin20_obj[["originalexp"]]@counts@x)) > 1e-8))
stopifnot(all(Matrix::colSums(bin20_obj[["originalexp"]]@counts) > 0))

# Create a clean RNA assay from originalexp counts
rna_counts <- GetAssayData(bin20_obj, assay = "originalexp", layer = "counts")
bin20_obj[["RNA"]] <- CreateAssayObject(counts = rna_counts)
DefaultAssay(bin20_obj) <- "RNA"


#Try to locate spatial coordinates in colData----
md <- bin20_obj@meta.data
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
                                   key = "SP_", assay = DefaultAssay(bin20_obj))
bin20_obj[["spatial"]] <- spatial_dr

# Basic scRNA-seq workflow----
plan(sequential)
bin20_obj <- NormalizeData(bin20_obj)
bin20_obj <- FindVariableFeatures(bin20_obj, nfeatures=2000)
bin20_obj <- ScaleData(bin20_obj, features = VariableFeatures(bin20_obj))

# 1) Make a sketch
set.seed(1)
n_ref <- 100000  # tune 100k–200k
ref_cells <- sample(colnames(bin20_obj), n_ref)
ref <- subset(bin20_obj, cells = ref_cells)

# 2) Do PCA/UMAP/clustering on the sketch
ref <- ScaleData(ref, features = VariableFeatures(bin20_obj), verbose = TRUE)
ref <- RunPCA(ref, features = VariableFeatures(bin20_obj),
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
  query = bin20_obj,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reduction = "pcaproject",    # <— key: project query onto ref PCs
  dims = 1:40,
  features = VariableFeatures(bin20_obj)
)

bin20_obj <- MapQuery(
  anchorset = anchors,
  query = bin20_obj,
  reference = ref,
  refdata = list(cluster = Idents(ref)),  # transfers cluster labels
  reference.reduction = "pca",
  reduction.model = "umap",               # projects into ref UMAP
  verbose = TRUE
)

# Results:
# - bin20_obj[["ref.pca"]]  : projected PCs for all cells
# - bin20_obj[["ref.umap"]] : projected UMAP for all cells
# - bin20_obj$predicted.cluster : transferred cluster labels

# 3B) If you prefer not to run anchors:
#     You can still project UMAP if you have PCs for all cells.
#     (e.g., if you later compute PCs, use:)
# emb_full <- uwot::umap_transform(
#   Embeddings(bin20_obj, "pca")[, 1:40],
#   model = ref@reductions$umap@misc$model
# )
# bin20_obj[["umap"]] <- CreateDimReducObject(embeddings = emb_full, key = "UMAP_")

Idents(bin20_obj) <-(as.numeric(bin20_obj@meta.data$predicted.cluster))

#save the original (uncleaned) object
library(SeuratDisk)
SaveH5Seurat(bin20_obj, filename = "20250906 Rat thymus stereo-seq recount/bin20_obj.h5seurat", overwrite = TRUE)
bin20_obj <- LoadH5Seurat("20250906 Rat thymus stereo-seq recount/bin20_obj.h5seurat") #Error: Call assays must have either a 'counts' or 'data' slot, missing for RNA

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

bin20_obj <- clean_for_save(bin20_obj)

#Set the leiden clusters as clusters of the Seurat object----
# 1) Make sure it's a clean factor (avoid weird types)
bin20_obj$leiden <- as.character(bin20_obj$leiden)
bin20_obj$leiden <- ifelse(is.na(bin20_obj$leiden), "NA", bin20_obj$leiden)
bin20_obj$leiden <- factor(bin20_obj$leiden)

# 2) Set active identities to SAW's Leiden
Idents(bin20_obj) <- "leiden"

bin20_obj$leiden <- factor(bin20_obj$leiden,
                             levels = sort(as.numeric(unique(bin20_obj$leiden))))

bin20_obj <- SetIdent(bin20_obj, value = factor(bin20_obj$leiden), levels = levels(unique(bin20_obj$leiden)))

# 3) Quick overview
table(Idents(bin20_obj))[1:20]   # first 20 clusters + sizes
length(levels(Idents(bin20_obj)))  # number of Leiden clusters

DimPlot(bin20_obj, reduction = "spatial", group.by = "leiden", label = TRUE)

# Save----
options(expressions = 5e5)               # raise recursion limit as a safety net
saveRDS(bin20_obj, "20250906 Rat thymus stereo-seq recount/bin20_obj.rds") 
bin20_obj <- readRDS("20250906 Rat thymus stereo-seq recount/bin20_obj.rds")

# Automatic mECA/mEFA recognition----
# pick signatures (adapt/expand as you like)
mTEC_genes  <- c("Epcam","Krt5","Krt8","Reg3a","Reg3g","Alox12e", "Fezf2")
mThymo_genes <- c("Cd69", "Tox", "Cd2", "Cd28", "Ms4a4a","Cd53","Sell","Ccr7")

# (1) module scores
bin20_obj <- AddModuleScore(bin20_obj, features=list(mTEC_genes),  name="mTEC_")
bin20_obj <- AddModuleScore(bin20_obj, features=list(mThymo_genes), name="mThymo_")

# (2) normalize scores (z within object)
md <- bin20_obj@meta.data
md$mTEC_norm  <- as.numeric(scale(md$mTEC_1))
md$mThymo_norm <- as.numeric(scale(md$mThymo_1))

# (3) simple rule to start:
#   mEFA:   mThymo high & mTEC low
#   mECA:   mThymo high & mTEC high
#   cortex:  everything else
q <- function(x,p=.75) as.numeric(quantile(x,p,na.rm=TRUE))

hi_mThymo <- md$mThymo_norm >= quantile(md$mThymo_norm,.90, na.rm=TRUE)
hi_mTEC   <- md$mTEC_norm  >= quantile(md$mTEC_norm,.85, na.rm=TRUE)
lo_mTEC   <- md$
  
  mTEC_norm  <= quantile(md$mTEC_norm,.84, na.rm=TRUE)

label <- ifelse(hi_mThymo & lo_mTEC, "mEFA",
                ifelse(hi_mThymo & hi_mTEC, "mECA", "cortex"))
bin20_obj$mRegion <- factor(label, levels=c("mEFA","mECA","cortex"))

a <-DimPlot(bin20_obj, reduction="spatial", group.by="mRegion", pt.size=1,
            cells     = rownames(bin20_obj@meta.data)[bin20_obj$mRegion %in% c("mEFA","mECA")]) +
  scale_y_reverse() + coord_fixed()

# (4) optional spatial smoothing (reduces salt-pepper):
# use 2D kNN average of scores
library(FNN)
k <- 200
knn <- get.knn(as.matrix(md[,c("x","y")]), k=k)$nn.index
smooth <- function(v) rowMeans(matrix(v[knn], ncol=k))
md$Thymo_smooth <- smooth(md$mThymo_norm); md$mTEC_smooth <- smooth(md$mTEC_norm)
hi_mThymo   <- md$Thymo_smooth >= q(md$Thymo_smooth,.90)
hi_mTEC   <- md$mTEC_smooth  >= q(md$mTEC_smooth,.90)
lo_mTEC   <- md$mTEC_smooth  <= q(md$mTEC_smooth,.89)
bin20_obj$mRegion <- factor(ifelse(hi_mThymo & lo_mTEC ,"mEFA", ifelse(hi_mThymo  & hi_mTEC ,"mECA","cortex")),
                              levels=c("mEFA","mECA","cortex"))

b <- DimPlot(
  bin20_obj,
  reduction = "spatial",
  group.by  = "mRegion",
  cells     = rownames(bin20_obj@meta.data)[bin20_obj$mRegion %in% c("mEFA","mECA")],
  pt.size   = 1
) +
  scale_y_reverse() + coord_fixed()

c <- DimPlot(
  bin20_obj,
  reduction = "spatial",
  group.by  = "mRegion",
  cells     = rownames(bin20_obj@meta.data)[bin20_obj$mRegion %in% c("mEFA","mECA")],
  pt.size   = 1
) +
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(15000, 11500))

a|b|c

DimPlot(bin20_obj, reduction = "spatial", group.by  = "mRegion",
  cols = c(mEFA="green", mECA="red",cortex="grey")
) + scale_y_reverse() + coord_fixed()
ggsave(filename = "pic/fig1C.pdf", plot = get_last_plot(), width = 7.4, height = 6.2, units = "in")

DimPlot(bin20_obj, reduction = "spatial", group.by  = "mRegion",
        cols = c(mEFA="green", mECA="red",cortex="grey"), pt.size = 4
        ) + scale_y_reverse() + coord_fixed(xlim = c(9000, 10300), ylim = c(7600, 4800)) + labs(title = NULL, x = NULL, y = NULL) + theme_void(base_size = 12)+ NoLegend() 
ggsave(filename = "pic/Fig1C_2.pdf", plot = get_last_plot(), width = 5.45, height = 6.28, units = "in")

DimPlot(bin20_obj, reduction = "spatial", group.by  = "mRegion",
        cols = c(mEFA="green", mECA="red",cortex="grey"), pt.size = 4
        ) + scale_y_reverse() + coord_fixed(xlim = c(11350, 12650), ylim = c(14600, 11800)) + labs(title = NULL, x = NULL, y = NULL) + theme_void(base_size = 12)+ NoLegend() 
ggsave(filename = "pic/Fig1C_3.pdf", plot = get_last_plot(), width = 5.45, height = 6.28, units = "in")

DimPlot(bin20_obj, reduction = "spatial", group.by  = "mRegion",
        cols = c(mEFA="green", mECA="red",cortex="grey"), pt.size = 4
        ) + scale_y_reverse() + coord_fixed(xlim = c(17300, 18600), ylim = c(6600, 3800)) + labs(title = NULL, x = NULL, y = NULL) + theme_void(base_size = 12)+ NoLegend() 
ggsave(filename = "pic/Fig1C_4.pdf", plot = get_last_plot(), width = 5.45, height = 6.28, units = "in")

# capsule recognition----
leiden13 <- subset(bin20_obj, subset = leiden == "13")

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

coords <- Embeddings(leiden13[["spatial"]])[,1:2]
df <- data.frame(cell = colnames(leiden13), x = coords[,1], y = coords[,2])

picked_cells <- lasso_select(df)
length(picked_cells)

leiden13$lasso_plotly <- FALSE
leiden13$lasso_plotly[match(picked_cells, colnames(leiden13))] <- TRUE

# Visualize the handwritten result
DimPlot(leiden13, reduction = "spatial")+
  scale_y_reverse() + coord_fixed() + labs(title = NULL) +  NoLegend() 
ggsave(filename = "pic/Fig1D_1.pdf", plot = get_last_plot(), width = 5.5, height = 5.5, units = "in")

DimPlot(leiden13, reduction = "spatial", group.by = "lasso_plotly")+
  scale_y_reverse() + coord_fixed() + labs(title = NULL) +  NoLegend() 
ggsave(filename = "pic/Fig1D_2.pdf", plot = get_last_plot(), width = 5.5, height = 5.5, units = "in")

# Save the handwritten result
saveRDS(leiden13, "20250813 Rat thymus stereo-seq analysis/capFb_bin20.rds")

leiden13 <- readRDS("20250813 Rat thymus stereo-seq analysis/capFb_bin20.rds")

# Map the lasso-selected capFbs back to the original object
picked_cells <- rownames(leiden13@meta.data[leiden13@meta.data$lasso_plotly == TRUE,])
bin20_obj$capsule_lasso <- FALSE
bin20_obj$capsule_lasso[match(picked_cells, colnames(bin20_obj))] <- TRUE

DimPlot(bin20_obj, reduction = "spatial", group.by = "capsule_lasso") +
  scale_y_reverse() + coord_fixed()

# Expand the selection to proximal bins
## 0) pre: bin20_obj$capsule_lasso is TRUE for your lasso picks
coords_all <- Embeddings(bin20_obj[["spatial"]])[, 1:2]
sel_idx    <- which(bin20_obj$capsule_lasso)

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

bin20_obj$capsule_expanded <- d_nearest <= radius

## 3) visualize
DimPlot(bin20_obj, reduction = "spatial", group.by = "capsule_expanded") +
  scale_y_reverse() + coord_fixed()

## 4) mix with mRegion
if (is.factor(bin20_obj$mRegion)) {
  levels(bin20_obj$mRegion) <- union(levels(bin20_obj$mRegion), "capsular")
}

bin20_obj$thymic_Regions <- as.character(bin20_obj$mRegion)
bin20_obj$capsule_dmin <- d_nearest

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

bin20_obj <- set_capsule_by_radius(bin20_obj, radius = 100)

bin20_obj$thyRegion <- factor(as.character(bin20_obj$thyRegion), levels = c("capsular", "cortex", "mECA", "mEFA"))

# Fig. 1D, 7.4 x 6.2 in
DimPlot(bin20_obj, reduction = "spatial", group.by = "thyRegion", 
        cols = c(capsular="royalblue", cortex="grey", mECA="red", mEFA="green"))+
  scale_y_reverse() + coord_fixed() + labs(title = NULL)
ggsave(filename = "pic/Fig1E.pdf", plot = get_last_plot(), width = 7.4, height = 6.2)

# Clustering and characterization 2 (RCTD)----
# https://satijalab.org/seurat/articles/spatial_vignette
RCTD_bin20 <- readRDS("20250906 Rat thymus stereo-seq recount/RCTD_bin20_20260114.rds")

# IDs RCTD actually used
rctd_ids <- rownames(RCTD_bin20@spatialRNA@coords)

# ---- build a compact metadata frame from your list-style results (top1/top2)
res_list <- RCTD_bin20@results

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
all_ids <- colnames(bin20_obj)
meta_all <- meta_rctd[intersect(rownames(meta_rctd), all_ids), , drop=FALSE]   # keep only valid columns
blank <- setdiff(all_ids, rownames(meta_all))
if (length(blank)) {
  na_block <- as.data.frame(matrix(NA, nrow=length(blank), ncol=ncol(meta_rctd)),
                            row.names = blank)
  names(na_block) <- names(meta_rctd)
  meta_all <- rbind(meta_all, na_block)
}
meta_all <- meta_all[all_ids, , drop=FALSE]  # order to match object

bin20_obj <- AddMetaData(bin20_obj, metadata = meta_all)

#Plots
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

bin20_obj@meta.data$RCTD_top1 <- factor(bin20_obj@meta.data$RCTD_top1, levels = lvl)
bin20_obj@meta.data$RCTD_top2 <- factor(bin20_obj@meta.data$RCTD_top2, levels = lvl)
bin20_obj@meta.data$RCTD_top3 <- factor(bin20_obj@meta.data$RCTD_top3, levels = lvl)

DimPlot(subset(bin20_obj, subset = RCTD_top2 %in% 
                 c("Projected: DN", "Projected: DP",  "Projected: Immature CD4", "Projected: Mature CD4",
                   "capFb", "TMC1", "TMC2", "TMC3", "TMC4")),
        reduction = "spatial", group.by = "RCTD_top2", split.by   = "RCTD_top2",
         ncol = 3) + 
  scale_y_reverse() + coord_fixed()+ NoLegend() + theme(plot.title = element_blank())
ggsave(filename = "pic/RCTD_result.pdf", plot = get_last_plot(), width = 8.6, height = 8.6)

DimPlot(subset(bin20_obj, subset = RCTD_top1_conf == T), reduction = "spatial",
  group.by = "RCTD_top1", split.by   = "RCTD_top1",
  pt.size    = 1, ncol = 7) + scale_y_reverse() + coord_fixed()

DimPlot(subset(bin20_obj, subset = RCTD_top2_conf == T), reduction = "spatial",
        group.by = "RCTD_top2", split.by   = "RCTD_top2",
        pt.size    = 1,ncol = 7) + scale_y_reverse() + coord_fixed()

DimPlot(subset(bin20_obj, subset = RCTD_top3_conf == T), reduction = "spatial",
        group.by = "RCTD_top3", split.by   = "RCTD_top3",
        pt.size    = 1, ncol = 7) + scale_y_reverse() + coord_fixed()

DimPlot(subset(bin20_obj, subset = RCTD_top1_conf == T), reduction = "spatial",
        group.by = "RCTD_top1", split.by   = "RCTD_top1",
        pt.size    = 3, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))

DimPlot(subset(bin20_obj, subset = RCTD_top2_conf == T), reduction = "spatial",
        group.by = "RCTD_top2", split.by   = "RCTD_top2",
        pt.size    = 3, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))

DimPlot(subset(bin20_obj, subset = RCTD_top3_conf == T), reduction = "spatial",
        group.by = "RCTD_top3", split.by   = "RCTD_top3",
        pt.size    = 2, ncol = 7) + 
  scale_y_reverse() + coord_fixed(xlim = c(11000, 13000), ylim = c(13500, 11500))

plot_rctd_bar <- function(obj_or_md, group_col="RCTD_top1", fill_col="leiden",
                          levels=NULL, K=10, label_frac=NULL,
                          legend_title=fill_col, y_limits=NULL) {
  library(dplyr); library(ggplot2); library(tibble)
  
  # get meta
  md <- if (inherits(obj_or_md, "Seurat")) obj_or_md@meta.data else obj_or_md
  md <- as_tibble(md, rownames = "cell")
  
  # tiny helper: force atomic character (handles list-columns)
  to_chr <- function(x) if (is.list(x))
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

plot_rctd_bar(subset(bin20_obj, subset = RCTD_top1_conf == T)@meta.data, group_col="RCTD_top1", levels=lvl, K=10, y_limits=c(0,80000))
plot_rctd_bar(subset(bin20_obj, subset = RCTD_top2_conf == T)@meta.data, group_col="RCTD_top2", levels=lvl, K=10, y_limits=c(0,80000))
plot_rctd_bar(subset(bin20_obj, subset = RCTD_top3_conf == T)@meta.data, group_col="RCTD_top3", levels=lvl, K=10, y_limits=c(0,80000))

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


plot_rctd_bar_simple(subset(bin20_obj, subset = RCTD_top1_conf == T)@meta.data, group_col="RCTD_top1", levels=lvl, y_limits=c(1,300000), y_log = T)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
ggsave(filename = "pic/figS3_5.pdf", plot = get_last_plot(), width = 7, height = 4.5, units = "in")

plot_rctd_bar_simple(subset(bin20_obj, subset = RCTD_top2_conf == T)@meta.data, group_col="RCTD_top2", levels=lvl, y_limits=c(1,300000), y_log = T)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 

plot_rctd_bar_simple(subset(bin20_obj, subset = RCTD_top3_conf == T)@meta.data, group_col="RCTD_top3", levels=lvl, y_limits=c(1,300000), y_log = T)+
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
plot_rctd_top_ratio(bin20_obj, top1_col = "RCTD_top1",  compare_col = "RCTD_top2",
                    ratio_label = "RCTD Top2/Top1")
ggsave(filename = "pic/figS3_6.pdf", plot = get_last_plot(), width = 7, height = 4.5, units = "in")

plot_rctd_top_ratio(bin20_obj, top1_col = "RCTD_top1",  compare_col = "RCTD_top3",
                    ratio_label = "RCTD Top3/Top1")
ggsave(filename = "pic/figS3_7.pdf", plot = get_last_plot(), width = 7, height = 4.5, units = "in")
# RCTD-predicted cell ratio in thymic regions----
plot_subset_ratio <- function(
    obj_or_md, target,
    layer_col   = "thymic_layers",
    top_cols    = c("RCTD_top1"),
    layer_lvls  = c("capsular","cortex","mECA","mEFA"),
    normalize   = c("by_layer"),        # or "by_target"
    conf        = TRUE,
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

plot_subset_ratio(bin20_obj, c("Projected: DN"))
plot_subset_ratio(bin20_obj, c("Projected: DP-P")) 
plot_subset_ratio(bin20_obj, c("Projected: DP")) 
plot_subset_ratio(bin20_obj, c("Projected: DP-Sig")) 
plot_subset_ratio(bin20_obj, c("Projected: Neg. Sel.")) 
plot_subset_ratio(bin20_obj, c("Projected: Immature CD4")) 
plot_subset_ratio(bin20_obj, c("Projected: Mature CD4")) 
plot_subset_ratio(bin20_obj, c("Projected: Immature CD8")) 
plot_subset_ratio(bin20_obj, c("Projected: Mature cycling T cell")) 
plot_subset_ratio(bin20_obj, c("capFb")) 
plot_subset_ratio(bin20_obj, c("TMC1")) 
plot_subset_ratio(bin20_obj, c("TMC2")) 
plot_subset_ratio(bin20_obj, c("TMC3")) 
plot_subset_ratio(bin20_obj, c("TMC4")) 
plot_subset_ratio(bin20_obj, c("vSMC_PC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:capEC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:vEC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:aEC")) 
plot_subset_ratio(bin20_obj, c("Projected: MEC"))
plot_subset_ratio(bin20_obj, c("Projected: TEC:early Pr")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:cTEC")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC1")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC-prol")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC2"))
plot_subset_ratio(bin20_obj, c("Projected: TEC:mimetic")) 

p1 <- plot_subset_ratio(bin20_obj, c("Projected: DN"))+theme(axis.title.x = element_blank())  
p2 <- plot_subset_ratio(bin20_obj, c("Projected: DP-P"))+theme(axis.title.x = element_blank())  
p3 <- plot_subset_ratio(bin20_obj, c("Projected: DP"))+theme(axis.title.x = element_blank())  
p4 <- plot_subset_ratio(bin20_obj, c("Projected: DP-Sig"))+theme(axis.title.x = element_blank())  
p5 <- plot_subset_ratio(bin20_obj, c("Projected: Neg. Sel."))+theme(axis.title.x = element_blank())  
p6 <- plot_subset_ratio(bin20_obj, c("Projected: Immature CD4"))+theme(axis.title.x = element_blank())  
p7 <- plot_subset_ratio(bin20_obj, c("Projected: Immature CD8"))+theme(axis.title.x = element_blank())  
p8 <- plot_subset_ratio(bin20_obj, c("Projected: Mature CD4"))+theme(axis.title.x = element_blank())
p9 <- plot_subset_ratio(bin20_obj, c("Projected: Treg"))+theme(axis.title.x = element_blank())  
p10 <- plot_subset_ratio(bin20_obj, c("Projected: Mature CD8"))+theme(axis.title.x = element_blank()) 

(p1|p2|p3|p4)/(p5|p6|p7)/(p8|p9|p10) 

plot_subset_ratio <- function(
    obj_or_md, target,
    layer_col   = "thymic_layers",
    top_cols    = c("RCTD_top1","RCTD_top2"),
    layer_lvls  = c("capsular","cortex","mECA","mEFA"),
    normalize   = c("by_layer"),        # or "by_target"
    conf        = TRUE,
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



plot_subset_ratio(bin20_obj, c("Projected: DN"))
plot_subset_ratio(bin20_obj, c("Projected: DP-P")) 
plot_subset_ratio(bin20_obj, c("Projected: DP")) 
plot_subset_ratio(bin20_obj, c("Projected: DP-Sig")) 
plot_subset_ratio(bin20_obj, c("Projected: Neg. Sel.")) 
plot_subset_ratio(bin20_obj, c("Projected: Immature CD4")) 
plot_subset_ratio(bin20_obj, c("Projected: Mature CD4")) 
plot_subset_ratio(bin20_obj, c("Projected: Immature CD8")) 
plot_subset_ratio(bin20_obj, c("Projected: Mature cycling T cell")) 
plot_subset_ratio(bin20_obj, c("capFb")) 
plot_subset_ratio(bin20_obj, c("TMC1")) 
plot_subset_ratio(bin20_obj, c("TMC2")) 
plot_subset_ratio(bin20_obj, c("TMC3")) 
plot_subset_ratio(bin20_obj, c("TMC4")) 
plot_subset_ratio(bin20_obj, c("vSMC_PC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:capEC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:vEC")) 
plot_subset_ratio(bin20_obj, c("Projected: EC:aEC")) 
plot_subset_ratio(bin20_obj, c("Projected: MEC"))
plot_subset_ratio(bin20_obj, c("Projected: TEC:early Pr")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:cTEC")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC1")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC-prol")) 
plot_subset_ratio(bin20_obj, c("Projected: TEC:mTEC2"))
plot_subset_ratio(bin20_obj, c("Projected: TEC:mimetic")) 

# Fig. 3B, 10 x 6 in
p1 <- plot_subset_ratio(bin20_obj, c("capFb"))+theme(axis.title.x = element_blank())
p2 <- plot_subset_ratio(bin20_obj, c("TMC1"))+theme(axis.title.x = element_blank())    
p3 <- plot_subset_ratio(bin20_obj, c("TMC2"))+theme(axis.title.x = element_blank())
p4 <- plot_subset_ratio(bin20_obj, c("TMC3"))+theme(axis.title.x = element_blank())
p5 <- plot_subset_ratio(bin20_obj, c("TMC4"))+theme(axis.title.x = element_blank())      
p6 <- plot_subset_ratio(bin20_obj, c("vSMC_PC"))

(p1|p2|p3)/(p4|p5|p6) 
ggsave(filename = "pic/Fig3B.pdf", plot = get_last_plot(), width = 10, height = 6)

p7 <- plot_subset_ratio(bin20_obj, c("Projected: EC:capEC")) 
p8 <- plot_subset_ratio(bin20_obj, c("Projected: EC:vEC")) 
p9 <- plot_subset_ratio(bin20_obj, c("Projected: EC:aEC")) 
(p7|p8|p9) 
ggsave(filename = "pic/Fig4B.pdf", plot = get_last_plot(), width = 10, height = 3.2)

DimPlot(subset(bin20_obj, subset=RCTD_top1 %in% c("capFb", "TMC1", "TMC2", "TMC3", "TMC4")), 
        reduction = "spatial", group.by = c("RCTD_top1"), pt.size = 2)+scale_y_reverse()|
  DimPlot(subset(bin20_obj, subset=RCTD_top2 %in% c("capFb", "TMC1", "TMC2", "TMC3", "TMC4")), 
          reduction = "spatial", group.by = c("RCTD_top2"), pt.size = 2)+scale_y_reverse() 

# gene expression comparison between mECA vs mEFA vs cortex vs capsule----
## read the reference dataset
RatThymicRef <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicRef.rds")
SimpleRatThymicRef.markers <- FindAllMarkers(JoinLayers(SimpleRatThymicRef), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SimpleRatThymicRef.markers <- SimpleRatThymicRef.markers %>% group_by(cluster) %>% slice_max(n = 100, order_by = avg_log2FC)

## Examine Fibroblast genes
DotPlot(SimpleRatThymicRef,  features = c("Pi16", "Cd34", "Pdgfra", "Pdgfrb", 
                                          "Crabp1", "Myoc", "Dpep1", "Obp3", "LOC102555263", # TMC1-2
                                          "Tshr", "Cela3b", "Gria4", "Apod", "Fam180a",      # TMC1-2
                                          "Alk", "Nrxn3", "Csmd1", "Nphs1", "LOC120103185",  # TMC3-4
                                           "LOC120095588", "Cdhr4", "Syt13", "Hdac9", "LOC120096691",
                                          "Pmfbp1", "Hrc", "Trpc6", "Npy1r", "Itih3",        # Pericyte
                                           "Sncg", "Fhl5", "Pcp4l1", "Scn3a", "Myh11"
                                    ))+RotatedAxis()

bin20_obj$group_fib <- factor(ifelse(bin20_obj$RCTD_top1 %in% c("capFb","TMC1","TMC2", "TMC3","TMC4", "vSMC_PC"), as.character(bin20_obj$RCTD_top1),
                                     ifelse(bin20_obj$RCTD_top2 %in% c("capFb","TMC1","TMC2", "TMC3","TMC4", "vSMC_PC"), as.character(bin20_obj$RCTD_top2), NA_character_)),
                              levels = c("capFb","TMC1","TMC2", "TMC3","TMC4", "vSMC_PC"))
fib_obj <- subset(bin20_obj, subset = !is.na(group_fib))

DotPlot(fib_obj, group.by = "group_fib", split.by = "thyRegion", cols = c("green", "blue", "red", "purple") , scale.max = 1,
        features = c("Pi16", "Cd34", "Pdgfra", "Pdgfrb", 
                       "Crabp1", "Myoc", "Dpep1", "Obp3", "LOC102555263",         # TMC1-2
                       "Tshr", "Cela3b", "Gria4", "Apod", "Fam180a",              # TMC1-2
                       "Alk", "Nrxn3", "Csmd1", "Nphs1", "LOC120103185",          # TMC3-4
                       "LOC120095588", "Cdhr4", "Syt13", "Hdac9", "LOC120096691", # TMC3-4
                       "Pmfbp1", "Hrc", "Trpc6", "Npy1r", "Itih3",                # Pericyte
                       "Sncg", "Fhl5", "Pcp4l1", "Scn3a", "Myh11")
        ) +RotatedAxis()+ 
  scale_y_discrete(limits = unlist(lapply(c("capFb","TMC1","TMC2", "TMC3","TMC4", "vSMC_PC"), function(g) paste(g, c("capsular","cortex","mECA","mEFA"), sep = "_"))))

DotPlot(bin20_obj,
        group.by = "thyRegion",
        features =  rev(c("Sphk1","Sphk2", "Spns2", "Sgpl1", "Plpp3"
        )), scale.max = NA)+
  RotatedAxis()+theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      legend.title = element_text(size = 10),
                      legend.text  = element_text(size = 9),
                      legend.key.height = unit(3, "mm"),
                      legend.key.width  = unit(3, "mm"))+
  coord_flip() +
  scale_size_area(
    max_size = 5,                         # visual max dot (adjust)
    limits   = c(0.01, 10),                 # clamp range
    oob      = squish,
    breaks   = c(0.01,0.1, 1, 10), labels = c("0.01", "0.1", "1", "10"))
#ggsave(filename = "pic/Fig3H_2.pdf", plot = get_last_plot(), width = 4.2, height = 3.2)

## Examine TEC genes
DotPlot(RatThymicRef,  features = c("Epcam", "Krt5", "Krt8",
                                    "Gpha2", "Anxa8", "Dapl1", "Krt13", "Dlk2", 
                                    "Epha5", "Serpina11", "LOC102548648", "Cntnap4", "LOC103691239",
                                    "Irx4", "Myo3a", "Ascl1", "Tspan33", 
                                    "Fezf2", "Calcb", "Skint1", "Reg3a", "Alox12e", "Wnt10a",
                                    "Tmprss13", "Gprc6a", "Bsnd", "Atp6v1g3", 
                                    "Wfdc10a", "Foxi1", "A2ml1", "Lgals5", "Asgr1", "Ppp1r17"
                                    ),
        group.by = "consolidated_cluster_s")+RotatedAxis()

bin20_obj$group_TEC <- factor(ifelse(bin20_obj$RCTD_top1 %in% c("Projected: TEC:early Pr", 
                                                                "Projected: TEC:cTEC", "Projected: TEC:mTEC1", "Projected: TEC:mTEC-prol", 
                                                                "Projected: TEC:mTEC2", "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic"), as.character(bin20_obj$RCTD_top1),
                                     ifelse(bin20_obj$RCTD_top2 %in% c("Projected: TEC:early Pr", 
                                                                       "Projected: TEC:cTEC", "Projected: TEC:mTEC1", "Projected: TEC:mTEC-prol", 
                                                                       "Projected: TEC:mTEC2", "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic"), as.character(bin20_obj$RCTD_top2), NA_character_)),
                              levels = c("Projected: TEC:early Pr", 
                                         "Projected: TEC:cTEC", "Projected: TEC:mTEC1", "Projected: TEC:mTEC-prol", 
                                         "Projected: TEC:mTEC2", "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic"))
TEC_obj <- subset(bin20_obj, subset = !is.na(group_TEC))

DotPlot(bin20_obj, group.by = "group_TEC", features = c("Epcam", "Krt5", "Krt8",
                                                       "Gpha2", "Anxa8", "Dapl1", "Krt13", "Dlk2", 
                                                       "Epha5", "Serpina11", "LOC102548648", "Cntnap4", "LOC103691239",
                                                       "Irx4", "Myo3a", "Ascl1", "Tspan33", 
                                                       "Fezf2", "Calcb", "Skint1", "Reg3a", "Alox12e", "Wnt10a",
                                                       "Tmprss13", "Gprc6a", "Bsnd", "Atp6v1g3", 
                                                       "Wfdc10a", "Foxi1", "A2ml1", "Lgals5", "Asgr1", "Ppp1r17"), scale.max = 1)+RotatedAxis()

## macrophage genes
DotPlot(RatThymicRef, scale.max = 2, features = c("Fcgr1a", "Itgam", "S100a8", "S100a9", "Cd163", "Flt3", "Itgax", "Zbtb46", "Xcr1", "Clec9a", "Cx3cr1", "Il1a"))+
  RotatedAxis()+theme(axis.title.x = element_blank())
FeaturePlot(RatThymicRef,  features = c("Fcgr1a", "Itgam", "S100a8", "S100a9", "Cd163", "Flt3", "Itgax", "Zbtb46", "Xcr1", "Clec9a", "Cx3cr1", "Il1a"))+
  RotatedAxis()+theme(axis.title.x = element_blank())

DotPlot(bin20_obj, group.by = "thyRegion", features = rev(c("Fcgr1a", "Itgam", "S100a8", "S100a9", "Cd163", "Flt3", "Itgax", "Zbtb46", "Xcr1", "Clec9a", "Cx3cr1", "Il1a")),
        scale.max = NA) + theme(axis.title = element_blank(), 
                                legend.title = element_text(size = 10),
                                legend.text  = element_text(size = 9),
                                legend.key.height = unit(3, "mm"),
                                legend.key.width  = unit(3, "mm")
                                ) + coord_flip() +
  scale_size_area(
    max_size = 5,                         # visual max dot (adjust)
    limits   = c(0.01, 7),                 # clamp range
    breaks   = c(0.1, 1, 5), labels = c("0.1", "1", "5"), trans="log1p")
ggsave(filename = "pic/Fig4E.pdf", plot = get_last_plot(), width = 5.2, height = 2.6)

## endothelial genes
FeaturePlot(RatThymicRef, reduction = "umap.PCA", features = c("Pecam1", "Cdh5", "Tek",
                                                               "Selp", "Vcam1", "Icam1", #ETP entrance
                                                               "Bst1", # T cell egress, PMID: 34276703
                                                               "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                                                               "Plpp3", #keep S1p concentration low
                                                               "Ccl19", "Ccl21", "Cxcl12" #egress-associated chemokine
                                                               ))

DotPlot(RatThymicRef,
        group.by = "consolidated_cluster",
               features = c("Pecam1", "Cdh5", "Tek",
                            "Selp", "Vcam1", "Icam1", #ETP entrance
                            "Bst1", # T cell egress, PMID: 34276703
                            "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                            "Plpp3", #keep S1p concentration low
                            "Ccl19", "Ccl21", "Cxcl12" #egress-associated chemokine
                            ))+RotatedAxis()+theme(axis.title.x = element_blank())

bin20_obj$group_Endo <- factor(ifelse(bin20_obj$RCTD_top1 %in% c("Projected: EC:capEC", "Projected: EC:aEC", "Projected: EC:vEC"),
                                      as.character(bin20_obj$RCTD_top1),
                                      ifelse(bin20_obj$RCTD_top2 %in% c("Projected: EC:capEC", "Projected: EC:aEC", "Projected: EC:vEC"), as.character(bin20_obj$RCTD_top2), NA_character_)),
                               levels = c("Projected: EC:capEC", "Projected: EC:aEC", "Projected: EC:vEC"))

Endo_obj <- subset(bin20_obj, subset = !is.na(group_Endo))

DotPlot(subset(Endo_obj, subset=thyRegion%in% c("mECA", "mEFA")), scale.max = 2,
        group.by = "thyRegion",
        features = c("Pecam1", "Cdh5", "Tek",
                     "Selp", "Vcam1", "Icam1", #ETP entrance
                     "Bst1", # T cell egress, PMID: 34276703
                     "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                     "Plpp3", #keep S1p concentration low
                     "Ccl19", "Ccl21", "Cxcl12" #egress-associated chemokine
        ))+RotatedAxis()+theme(axis.title.x = element_blank())

# Apply mECA/mEFA-EC score to EC subsets in reference scRNA-seq
agg_means <- AverageExpression(Endo_obj , group.by="mRegion", assays="RNA",  slot="data")

ranks <- (agg_means$RNA[,"mECA"])-(agg_means$RNA[,"mEFA"])
names(ranks) <- rownames(agg_means$RNA)
ranks <- sort(ranks, decreasing=TRUE)

# make reference marker table for subtraction
RatThymicRef <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicRef.rds")
RatThymicRef <- UpdateSeuratObject(RatThymicRef)
SimpleRatThymicRef <- RenameIdents(RatThymicRef, "Early SP"="thymocyte", "Mature SP"="thymocyte",
                                   "TMC1"="TMC", "TMC2"="TMC", "TMC3"="TMC","TMC4"="TMC", "vSMC/PC"="TMC",
                                   "TEC-a"="cTEC", "TEC-b"="cTEC", "TEC-c"="cTEC", "TEC-d"="mTEC", "TEC-e"="mimetic cells")

DimPlot(SimpleRatThymicRef, reduction = "umap.PCA", label = T, label.box = T, repel = T) +NoLegend()

RatThymicMarker <- FindAllMarkers(JoinLayers(SimpleRatThymicRef), logfc.threshold = 1, min.pct = 0.1, only.pos = T)
RatThymicMarker <- RatThymicMarker %>% filter(avg_log2FC > 1, pct.1 > 0.1) %>% arrange(cluster, desc(avg_log2FC))

bad_genes <- unique(na.omit(RatThymicMarker$gene[ RatThymicMarker$cluster %in% c("mTEC", "TMC", "thymocyte")]))
keep <- !(names(ranks) %in% bad_genes | grepl("^(Rpl|Rps)", names(ranks)))
ranks_TMC_mTECdeleted <- ranks[keep]

# mECA thymocyte vs mEFA thymocyte scores in reference scRNA-seq
mECAEndoMarker <- names(ranks_TMC_mTECdeleted )[order(ranks_TMC_mTECdeleted , decreasing = TRUE)[seq_len(min(100, length(ranks_TMC_mTECdeleted )))]]
mEFAEndoMarker <- names(ranks_TMC_mTECdeleted )[order(ranks_TMC_mTECdeleted , decreasing = F)[seq_len(min(100, length(ranks_TMC_mTECdeleted )))]]

RatThymicEndo.obj <- subset(RatThymicStroma.obj, ident=c("Endo-1", "Endo-2"))
RatThymicEndo.obj <- RunUMAP (RatThymicEndo.obj, dims = 1:17)
RatThymicEndo.obj <- FindClusters(RatThymicEndo.obj, resolution = 0.4)

RatThymicEndo.obj <- AddModuleScore(
  object   = RatThymicEndo.obj,
  features = list(mECAEndoMarker),
  name     = "mECA Endothelial score")

RatThymicEndo.obj <- AddModuleScore(
  object   = RatThymicEndo.obj,
  features = list(mEFAEndoMarker),
  name     = "mEFA Endothelial score")

FeaturePlot(RatThymicEndo.obj, features = c("mECA Endothelial score1", "mEFA Endothelial score1"))

DotPlot(RatThymicEndo.obj, features = c("mECA Endothelial score1", "mEFA Endothelial score1"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(Endo_obj, group.by = "group_Endo", split.by = "thyRegion", cols = c("green", "blue", "red", "purple"), features = c("Pecam1", "Cdh5", "Tek",
                                                       "Selp", "Vcam1", "Icam1", #ETP entrance
                                                       "Bst1", # T cell egress
                                                       "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                                                       "Plpp3", #keep S1p concentration low
                                                       "Ccl19", "Ccl21", "Cxcl12" #egress-associated chemokine
                                                       ), scale.max = 10)+RotatedAxis()+theme(axis.title.x = element_blank())

DotPlot(Endo_obj, group.by = "thyRegion", features = RatThymicEndo.markers %>% 
          filter(cluster %in% 0:6) %>% group_by(cluster) %>%  slice_head(n = 20) %>% pull(gene),
        scale.max = 1)+RotatedAxis()+theme(axis.title.x = element_blank())


RatThymicEndo <- subset(RatThymicRef, idents ="Endo")

RatThymicEndo <- JoinLayers(RatThymicEndo)

RatThymicEndo <- FindClusters(RatThymicEndo, resolution = 0.4)
RatThymicEndo <- RunUMAP(RatThymicEndo, dims = 1:15, reduction = "pca", reduction.name = "Endo.umap")

DimPlot(RatThymicEndo, reduction = "Endo.umap", label = T, label.box = T)

DotPlot(RatThymicEndo, features = c("Pecam1", "Cdh5", "Selp", "Vcam1", "Icam1", #ETP entrance
                                    "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                                    "Plpp3", #keep S1p concentration low
                                    "Ccl19", "Ccl21", "Cxcl12"  #egress-associated chemokine
                                    ))+RotatedAxis()

DotPlot(bin20_obj,  group.by = "thyRegion", features = c("Pecam1", "Cdh5", "Selp", "Vcam1", "Icam1", #ETP entrance
                                    "Spns2", "Sphk1", "Sphk2",  #Mature T cell egress
                                    "Plpp3", #keep S1p concentration low
                                    "Ccl19", "Ccl21", "Cxcl12"  #egress-associated chemokine
                                ), scale.max = 0.1)+RotatedAxis()

RatThymicEndo.markers <- FindAllMarkers(RatThymicEndo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
RatThymicEndo.markers <- RatThymicEndo.markers %>% group_by(cluster) %>% slice_max(n = 100, order_by = avg_log2FC)

DotPlot(RatThymicRef, features = c("Pecam1", "Cdh5", #Total Endothelial cells
                                   "Slc26a10", "Nxph1", #Cluster0
                                   "Selp", "Sele", #Cluster1
                                   "Car4", "C1qtnf9" #Cluster2
                           ))+RotatedAxis()

DotPlot(bin20_obj, group.by = "thyRegion", features = c("Pecam1", "Cdh5", #Total Endothelial cells
                                   "Slc26a10", "Nxph1", #Cluster0
                                   "Selp", "Sele", #Cluster1
                                   "Car4", "C1qtnf9"#Cluster2
                                   ), scale.max = 0.5)+RotatedAxis()

# above all, ETP-entrance blood vessels are located in mECA (Selp^+).
# But localization of T cell-egress blood vessel is unknown, because the related genes (expecially BST-1) are expressed not only from endothelial cells but also from other cells.
# and EC bins are contaminated with other cell-derived mRNAs.