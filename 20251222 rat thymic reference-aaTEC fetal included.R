#https://satijalab.org/seurat/articles/pbmc3k_tutorial
library(Seurat)
library(dplyr)
library(ggplot2)
library(orthogene)
# Rat thymocyte object----
RatThymocyte.obj <- CreateSeuratObject(counts = Read10X(
  "20250904 Rat thymocyte scRNA-seq/filtered_feature_bc_matrix/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
), project = "RatThymocyte")
RatThymocyte.obj

RatThymocyte.obj[["percent.mt"]] <- PercentageFeatureSet(RatThymocyte.obj, pattern = "^mt-")
head(RatThymocyte.obj@meta.data, 5)
VlnPlot(RatThymocyte.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mean(RatThymocyte.obj@meta.data$nFeature_RNA)

plot1 <- FeatureScatter(RatThymocyte.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RatThymocyte.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

RatThymocyte.obj <- subset(RatThymocyte.obj, subset = nFeature_RNA > 200)

RatThymocyte.obj <- NormalizeData(RatThymocyte.obj)

RatThymocyte.obj <- FindVariableFeatures(RatThymocyte.obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RatThymocyte.obj), 10)

plot1 <- VariableFeaturePlot(RatThymocyte.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(RatThymocyte.obj)
RatThymocyte.obj <- ScaleData(RatThymocyte.obj, features = all.genes)


RatThymocyte.obj <- RunPCA(RatThymocyte.obj, features = VariableFeatures(object = RatThymocyte.obj))
print(RatThymocyte.obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(RatThymocyte.obj, dims = 1:2, reduction = "pca")
DimPlot(RatThymocyte.obj, reduction = "pca")

DimHeatmap(RatThymocyte.obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(RatThymocyte.obj, dims = 1:15, cells = 500, balanced = TRUE)

RatThymocyte.obj <- JackStraw(RatThymocyte.obj, num.replicate = 100)
RatThymocyte.obj <- ScoreJackStraw(RatThymocyte.obj, dims = 1:15)
JackStrawPlot(RatThymocyte.obj, dims = 1:15)
ElbowPlot(RatThymocyte.obj)

RatThymocyte.obj <- FindNeighbors(RatThymocyte.obj)
RatThymocyte.obj <- FindClusters(RatThymocyte.obj, resolution = 0.5)
head(Idents(RatThymocyte.obj), 5)

RatThymocyte.obj <- RunUMAP (RatThymocyte.obj, dims = 1:15)

DimPlot(RatThymocyte.obj, label = TRUE, label.box = T, repel = T, reduction = "umap")+NoLegend()

saveRDS(RatThymocyte.obj, file = "20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

# Rat thymocyte: Cluster identification
## PMID: 37580604
DotPlot(RatThymocyte.obj, features = c("Snca", "Alas2", "Hba-a3", "Rsad2", "Oasl2", #Erythrocyte
                                       "ENSRNOG00000064589", "Ms4a1", "Cd19", #B cells
                                       "Ctsh", "Spi1", "Tyrobp", "Mpeg1", #Myeloid
                                       "Cpa3", "Hes1", #GD T cells
                                       "Mpzl2",  #DN
                                       "C7h12orf75", "H2ac10", "Esco2", "Mki67", #DP (P), Mature cycling
                                       "Slc7a11", #DP (Q1)
                                       "P2rx1", "Rag1", #DP (Q2)
                                       "Rapsn", "Ccr4", #DP (Sig)
                                       "Ifit1", "RGD1309362", #Interferon sig.
                                       "Itm2a", "Zbtb7b", "Cd40lg", "Chl1", "Gpr83", "Acvrl1", #Immature CD4, Neg. sel(2)
                                       "Tesc",  #Immature CD8
                                       "S1pr1", "Slfn1", #Mature CD4
                                       "Cd7", #NKT
                                       "Dapl1",#Mature CD8
                                       "Nkg7", #NKT
                                       "Egr2", "Nr4a3", "Pdcd1lg2", #Neg. sel. (1)
                                       "Tnfrsf4", "Tnfrsf9", "Arhgap20", "Foxp3", "Slc22a2", "LOC120102904", #Treg
                                       "S100a4", "Zbtb16", "Klra22")  #NKT
)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(RatThymocyte.obj, features = c("Notch1", "Notch2", "Hes1", # PMID: 36316479, Notch signaling
                                       "Mpzl2", #PMID: 37580604, DN
                                       "H2ac10", "Esco2", "Mki67", #DP (P)
                                       "Slc7a11", "Dntt", "Rag1", "Rag2", "Thy1", "Olfml3", "Ada", # DP(P)~DP)
                                       "Ccdc172", "Gucy1b2", "Sox6","Pde6c", "Tecta", "Spdl1", # DP
                                       "Tox2", "Hivep3", "Itm2a", "Nab2", "Srgn", "Tox", "Nr4a1", "Cd28", "Rgs2", "Cd69", "Raph1", "Cd2", # Early SP
                                       "Klrk1", "Ms4a4c", "Calhm6", "Cd79a", "Sell", "Ccnd2", "Cd40lg", "Rasa3", "Ms4a4a", "Slc4a4", "S1pr1",
                                       "Cd96", "Ccr7", "Foxp3", "Il2ra" # Mature SP
))+RotatedAxis()

# Rat Thymocyte: Name clusters
RatThymoNames <- c("DP1", "DP2", "Early SP", "Mature SP", "DP-P1", "DP-P2", "DP1", "DP-P3", "B+myeloid", "DP2", "DN",  "B+myeloid")
names(RatThymoNames) <- levels(Idents(RatThymocyte.obj))
RatThymocyte.obj <- RenameIdents(object = RatThymocyte.obj, RatThymoNames)

DimPlot(RatThymocyte.obj, label = T, label.box = T, repel = T)+NoLegend()

# Define some subsets by handwriting
## Tiny Shiny “gadget”
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
          yaxis=list(title="SP_2")
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

coords <- Embeddings(RatThymocyte.obj[["umap"]])[,1:2]
df <- data.frame(cell = colnames(RatThymocyte.obj), x = coords[,1], y = coords[,2])

## add DN cells
picked_cells <- lasso_select(df)

RatThymocyte.obj$hand_cluster <- as.character(Idents(RatThymocyte.obj))
RatThymocyte.obj$hand_cluster[match(picked_cells, colnames(RatThymocyte.obj))] <- "DN"
Idents(RatThymocyte.obj) <- RatThymocyte.obj$hand_cluster

## add Mature SP cells
picked_cells <- lasso_select(df)

RatThymocyte.obj$hand_cluster[match(picked_cells, colnames(RatThymocyte.obj))] <- "Mature SP"
Idents(RatThymocyte.obj) <- RatThymocyte.obj$hand_cluster


RatThymo_order <- c("B+myeloid", "DN", "DP-P1", "DP-P2", "DP-P3", "DP1", "DP2", "Early SP", "Mature SP")

RatThymocyte.obj <- SetIdent(RatThymocyte.obj, value = factor(Idents(RatThymocyte.obj), levels = RatThymo_order))

DimPlot(RatThymocyte.obj, label = TRUE, label.box = T, repel = T)

saveRDS(RatThymocyte.obj, file = "20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

DimPlot(RatThymocyte.obj, label = TRUE, label.box = T, repel = T)
ggsave(filename = "pic/FigS1A.pdf", plot = get_last_plot(), width = 7, height = 5)

DotPlot(RatThymocyte.obj, features = c("Cd19", # B cells
                                       "Ctsh", "Spi1", "Tyrobp", # Myeloid
                                       "Notch1", "Hes1", # PMID: 36316479, Notch signaling
                                       "Mpzl2", #PMID: 37580604, DN
                                       "H2ac10", "Esco2", "Mki67", #DP (P)
                                       "Slc7a11", "Rag1", "Rag2", "Thy1", # DP(P)~DP
                                       "Ccdc172", "Gucy1b2", "Sox6", # DP
                                       "Tox", "Tox2", "Itm2a", "Nr4a1", "Cd2", "Cd28",  "Cd69",  # Early SP
                                       "Ms4a4a", "Sell", "Cd40lg", "Cd79a", "Cd96", "S1pr1", "Ccr7"# Mature SP
                                       ))+
  scale_y_discrete(limits = rev(levels(RatThymocyte.obj))) +
  RotatedAxis() + labs(x=NULL, y=NULL)
ggsave(filename = "pic/FigS1A_2.pdf", plot = get_last_plot(), width = 10, height = 4.2)

# Read mouse thymocytes and project labels to rat thymocytes----
ReadMouseThymocytes <- function(
    dir,
    sample_id,
    project        = sample_id,
    gene_assay_key = "Gene Expression",
    adt_assay_key  = "Antibody Capture",
    adt_name       = "ADT",
    gene.column    = 2,
    cell.column    = 1,
    unique.features = TRUE,
    strip.suffix    = FALSE
) {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  # We'll need SeuratObject for v5's CreateAssay
  requireNamespace("SeuratObject", quietly = TRUE)
  
  x <- Seurat::Read10X(
    data.dir        = dir,
    gene.column     = gene.column,
    cell.column     = cell.column,
    unique.features = unique.features,
    strip.suffix    = strip.suffix
  )
  
  if (!gene_assay_key %in% names(x)) {
    stop(sprintf("Could not find '%s' in Read10X() output. Available: %s",
                 gene_assay_key, paste(names(x), collapse = ", ")))
  }
  gene <- x[[gene_assay_key]]
  adt  <- if (!is.null(adt_assay_key) && adt_assay_key %in% names(x)) x[[adt_assay_key]] else NULL
  rownames(adt) <- rownames(adt) %>% gsub("^ADT_", "", .) %>% gsub("_A[0-9]+$", "", .) %>% gsub("_TotalA", "", .) %>% gsub("_", "", .)
  
  # Base object with RNA
  so <- Seurat::CreateSeuratObject(counts = gene, project = project)
  
  # Add ADT if present — handle Seurat v5 vs v4 gracefully
  if (!is.null(adt)) {
    if (utils::packageVersion("SeuratObject") >= "5.0.0" &&
        "CreateAssay" %in% getNamespaceExports("SeuratObject")) {
      # Seurat v5: create Assay5 via CreateAssay(layers = list(counts = ...))
      so[[adt_name]] <- SeuratObject::CreateAssay(layers = list(counts = adt), key = paste0(adt_name, "_"))
    } else {
      # Seurat v4 fallback
      so[[adt_name]] <- Seurat::CreateAssayObject(counts = adt)
    }
  }
  
  # Prefix cell IDs
  so <- Seurat::RenameCells(so, add.cell.id = sample_id)
  so
}

# B6 mouse
B6_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631537_ZRS06_1_B6", sample_id = "ZRS06_1_B6", project = "ZRS06_1_B6")
B6_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631540_ZRS07_1_B6_r2", sample_id = "ZRS07_1_B6_r2", project = "ZRS07_1_B6_r2")
B6_r3 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631544_ZRS07_5_B6_r3", sample_id = "ZRS07_5_B6_r3", project = "ZRS07_5_B6_r3")
B6_r4_lane1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631548_1_B6_thy_r4_111", sample_id = "1_B6_thy_r4_111", project = "1_B6_thy_r4_111")
B6_r4_lane2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631549_2_B6_thy_r4_111", sample_id = "2_B6_thy_r4_111", project = "2_B6_thy_r4_111")
B6_r4 <- merge(B6_r4_lane1, B6_r4_lane2, project ="B6_thy_r4_111")
B6_r5_lane8 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631552_8_B6_thy_r5_111", sample_id = "8_B6_thy_r5_111", project = "8_B6_thy_r5_111")
B6_r5_lane9 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631553_9_B6_thy_r5_111", sample_id = "9_B6_thy_r5_111", project = "9_B6_thy_r5_111")
B6_r5 <- merge(B6_r5_lane8, B6_r5_lane9, project ="B6_thy_r8_111")

B6_thymocyte <- merge(B6_r1, c(B6_r2, B6_r3, B6_r4, B6_r5), project ="mouse_thymocyte")

B6_thymocyte[["percent.mt"]] <- PercentageFeatureSet(B6_thymocyte, pattern = "^mt-")
head(B6_thymocyte@meta.data, 5)
VlnPlot(B6_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

B6_thymocyte <- subset(B6_thymocyte, subset = nFeature_RNA > 200)

B6_thymocyte@meta.data <- B6_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("ZRS06_1_B6", orig.ident) ~ "B6_r1",
    grepl("ZRS07_1_B6_r2", orig.ident) ~ "B6_r2",
    grepl("ZRS07_5_B6_r3", orig.ident) ~ "B6_r3",
    grepl("^1_B6_thy_r4_111", orig.ident) | grepl("^2_B6_thy_r4_111", orig.ident) ~ "B6_r4",
    grepl("^8_B6_thy_r5_111", orig.ident) | grepl("^9_B6_thy_r5_111", orig.ident) ~ "B6_r5",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#AND mouse
AND_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631538_ZRS06_2_AND", sample_id = "ZRS06_2_AND", project = "ZRS06_2_AND")
AND_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631543_ZRS07_4_AND_r2", sample_id = "ZRS07_4_AND_r2", project = "ZRS07_4_AND_r2")

AND_thymocyte <- merge(AND_r1, AND_r2, project ="AND_thymocyte")

AND_thymocyte[["percent.mt"]] <- PercentageFeatureSet(AND_thymocyte, pattern = "^mt-")
head(AND_thymocyte@meta.data, 5)
VlnPlot(AND_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

AND_thymocyte <- subset(AND_thymocyte, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 9)

AND_thymocyte@meta.data <- AND_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("ZRS06_2_AND", orig.ident) ~ "AND_r1",
    grepl("ZRS07_4_AND_r2", orig.ident) ~ "AND_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#F5 mouse
F5_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631539_ZRS06_3_F5", sample_id = "ZRS06_3_F5", project = "ZRS06_3_F5")


F5_thymocyte <- merge(F5_r1, F5_r2, project ="F5_thymocyte")

F5_thymocyte[["percent.mt"]] <- PercentageFeatureSet(F5_thymocyte, pattern = "^mt-")
head(F5_thymocyte@meta.data, 5)
VlnPlot(F5_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

F5_thymocyte <- subset(F5_thymocyte, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 7)

F5_thymocyte@meta.data <- F5_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("ZRS06_3_F5", orig.ident) ~ "F5_r1",
    grepl("ZRS07_8_F5_r2", orig.ident) ~ "F5_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#B2M (MHCI_KO) mouse
MHCI_KO_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631541_ZRS07_2_B2M_r1", sample_id = "ZRS07_2_B2M_r1", project = "ZRS07_2_B2M_r1")
MHCI_KO_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631542_ZRS07_3_B2M_r2", sample_id = "ZRS07_3_B2M_r2", project = "ZRS07_3_B2M_r2")

MHCI_KO_thymocyte <- merge(MHCI_KO_r1, MHCI_KO_r2, project ="MHCI_KO_thymocyte")

MHCI_KO_thymocyte[["percent.mt"]] <- PercentageFeatureSet(MHCI_KO_thymocyte, pattern = "^mt-")
head(MHCI_KO_thymocyte@meta.data, 5)
VlnPlot(MHCI_KO_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MHCI_KO_thymocyte <- subset(MHCI_KO_thymocyte, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 9)

MHCI_KO_thymocyte@meta.data <- MHCI_KO_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("ZRS07_2_B2M_r1", orig.ident) ~ "MHCI_KO_r1",
    grepl("ZRS07_3_B2M_r2", orig.ident) ~ "MHCI_KO_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#MHCII_KO mouse
MHCII_KO_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631545_ZRS07_6_MHC2_r1", sample_id = "ZRS07_6_MHC2_r1", project = "ZRS07_6_MHC2_r1")
MHCII_KO_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631546_ZRS07_7_MHC2_r2", sample_id = "ZRS07_7_MHC2_r2", project = "ZRS07_7_MHC2_r2")

MHCII_KO_thymocyte <- merge(MHCII_KO_r1, MHCII_KO_r2, project ="MHCII_KO_thymocyte")

MHCII_KO_thymocyte[["percent.mt"]] <- PercentageFeatureSet(MHCII_KO_thymocyte, pattern = "^mt-")
head(MHCII_KO_thymocyte@meta.data, 5)
VlnPlot(MHCII_KO_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MHCII_KO_thymocyte <- subset(MHCII_KO_thymocyte, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 9)

MHCII_KO_thymocyte@meta.data <- MHCII_KO_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("ZRS07_6_MHC2_r1", orig.ident) ~ "MHCII_KO_r1",
    grepl("ZRS07_7_MHC2_r2", orig.ident) ~ "MHCII_KO_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#OT_I mouse
OT_I_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631550_3_OT1_thy_r1_111", sample_id = "3_OT1_thy_r1_111", project = "3_OT1_thy_r1_111")
OT_I_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631554_10_OT1_thy_r2_111", sample_id = "10_OT1_thy_r2_111", project = "10_OT1_thy_r2_111")

OT_I_thymocyte <- merge(OT_I_r1, OT_I_r2, project ="OT_I_thymocyte")

OT_I_thymocyte[["percent.mt"]] <- PercentageFeatureSet(OT_I_thymocyte, pattern = "^mt-")
head(OT_I_thymocyte@meta.data, 5)
VlnPlot(OT_I_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

OT_I_thymocyte <- subset(
  OT_I_thymocyte,
  subset =
    (orig.ident == "3_OT1_thy_r1_111"  & nFeature_RNA > 1500 & nFeature_RNA < 5000 & percent.mt < 9) |
    (orig.ident == "10_OT1_thy_r2_111" & nFeature_RNA > 700  & nFeature_RNA < 3500 & percent.mt < 9)
)

OT_I_thymocyte@meta.data <- OT_I_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("3_OT1_thy_r1_111", orig.ident) ~ "OT_I_r1",
    grepl("10_OT1_thy_r2_111", orig.ident) ~ "OT_I_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

#OT_II mouse
OT_II_r1 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631551_4_OT2_thy_r1_111", sample_id = "4_OT2_thy_r1_111", project = "4_OT2_thy_r1_111")
OT_II_r2 <- ReadMouseThymocytes(dir = "PMID_37580604 mouse thymocyte/GSM5631555_11_OT2_thy_r2_111", sample_id = "11_OT2_thy_r2_111", project = "11_OT2_thy_r2_111")

OT_II_thymocyte <- merge(OT_II_r1, OT_II_r2, project ="OT_II_thymocyte")

OT_II_thymocyte[["percent.mt"]] <- PercentageFeatureSet(OT_II_thymocyte, pattern = "^mt-")
head(OT_II_thymocyte@meta.data, 5)
VlnPlot(OT_II_thymocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

OT_II_thymocyte <- subset(
  OT_II_thymocyte,
  subset =
    (orig.ident == "4_OT2_thy_r1_111"  & nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 9) |
    (orig.ident == "11_OT2_thy_r2_111" & nFeature_RNA > 500  & nFeature_RNA < 2500 & percent.mt < 9)
)

OT_II_thymocyte@meta.data <- OT_II_thymocyte@meta.data %>%
  mutate(individual = case_when(
    grepl("4_OT2_thy_r1_111", orig.ident) ~ "OT_II_r1",
    grepl("11_OT2_thy_r2_111", orig.ident) ~ "OT_II_r2",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))


#merge all genotypes
all_m_thymocyte <- merge(B6_thymocyte, 
                         c(AND_thymocyte, F5_thymocyte, MHCI_KO_thymocyte, MHCII_KO_thymocyte, OT_I_thymocyte, OT_II_thymocyte),
                         project ="all_mouse_thymocyte")

all_m_thymocyte@meta.data <- all_m_thymocyte@meta.data %>%
  mutate(genotype = case_when(
    grepl("B6", individual) ~ "WT",
    grepl("AND", individual) ~ "AND",
    grepl("F5", individual) ~ "F5",
    grepl("MHCI_KO", individual) ~ "MHCI_KO",
    grepl("MHCII_KO", individual) ~ "MHCII_KO",
    grepl("OT_I_r", individual) ~ "OT_I",
    grepl("OT_II", individual) ~ "OT_II",
    TRUE ~ NA_character_  # fallback for unmatched entries
  ))

all_m_thymocyte

all_m_thymocyte <- all_m_thymocyte %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 2, cluster.name = "unintegrated_clusters")
all_m_thymocyte <- NormalizeData(all_m_thymocyte, normalization.method = "CLR", margin = 2, assay = "ADT")

all_m_thymocyte <- RunUMAP(all_m_thymocyte, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(all_m_thymocyte, reduction = "umap.unintegrated", split.by = "genotype")

saveRDS(all_m_thymocyte, "PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")

all_m_thymocyte <- IntegrateLayers(all_m_thymocyte, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
all_m_thymocyte <- FindNeighbors(all_m_thymocyte, reduction = "integrated.cca", dims = 1:30, assay = "RNA")
all_m_thymocyte <- FindClusters(all_m_thymocyte, resolution = 1, cluster.name = "cca_clusters", assay = "RNA")
all_m_thymocyte  <- RunUMAP(all_m_thymocyte , reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(all_m_thymocyte, reduction = "umap.cca", split.by = "genotype")
DimPlot(all_m_thymocyte, reduction = "umap.cca", label = T, label.box = T)

all_m_thymocyte <- IntegrateLayers(all_m_thymocyte, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", assay = "RNA")
all_m_thymocyte <- FindNeighbors(all_m_thymocyte, reduction = "integrated.rpca", dims = 1:30, assay = "RNA")
all_m_thymocyte <- FindClusters(all_m_thymocyte, resolution = 1, cluster.name = "rpca_clusters", assay = "RNA")
all_m_thymocyte  <- RunUMAP(all_m_thymocyte , reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", assay = "RNA")
DimPlot(all_m_thymocyte, reduction = "umap.rpca", split.by = "genotype")
DimPlot(all_m_thymocyte, reduction = "umap.rpca", label = T, label.box = T)

all_m_thymocyte <- IntegrateLayers(all_m_thymocyte, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
all_m_thymocyte <- FindNeighbors(all_m_thymocyte, reduction = "harmony", dims = 1:30, assay = "RNA")
all_m_thymocyte <- FindClusters(all_m_thymocyte, resolution = 1, cluster.name = "harmony_clusters", assay = "RNA")
all_m_thymocyte  <- RunUMAP(all_m_thymocyte , reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(all_m_thymocyte, reduction = "umap.harmony", split.by = "genotype")
DimPlot(all_m_thymocyte, reduction = "umap.harmony", label = T, label.box = T)

all_m_thymocyte <- IntegrateLayers(all_m_thymocyte, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.PCA")
all_m_thymocyte <- FindNeighbors(all_m_thymocyte, reduction = "integrated.PCA", dims = 1:30, assay = "RNA")
all_m_thymocyte <- FindClusters(all_m_thymocyte, resolution = 1, cluster.name = "PCA_clusters", assay = "RNA")
all_m_thymocyte  <- RunUMAP(all_m_thymocyte , reduction = "integrated.PCA", dims = 1:30, reduction.name = "umap.PCA")
DimPlot(all_m_thymocyte, reduction = "umap.PCA", split.by = "genotype")
DimPlot(all_m_thymocyte, reduction = "umap.PCA", label = T, label.box = T)

saveRDS(all_m_thymocyte, "PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")
all_m_thymocyte <- readRDS("PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")

#check integration results
DimPlot(all_m_thymocyte, reduction = "umap.cca", split.by = "genotype", group.by = "cca_clusters")
FeaturePlot(all_m_thymocyte, 
            features = c("Cd24a", "Cd69", "Sell", "Cd4", "Cd8b1", "Foxp3", "Tcrg-C2", "Mki67", "Mpzl2", "1500009L16Rik"),
            reduction = "umap.cca")

DimPlot(all_m_thymocyte, reduction = "umap.rpca", split.by = "genotype", group.by = "rpca_clusters")
FeaturePlot(all_m_thymocyte, 
            features = c("Cd24a", "Cd69", "Sell", "Cd4", "Cd8b1", "Foxp3", "Tcrg-C2", "Mki67", "Mpzl2", "1500009L16Rik"),
            reduction = "umap.rpca")

DimPlot(all_m_thymocyte, reduction = "umap.harmony", split.by = "genotype", group.by = "harmony_clusters")
FeaturePlot(all_m_thymocyte, 
            features = c("Cd24a", "Cd69", "Sell", "Cd4", "Cd8b1", "Foxp3", "Tcrg-C2", "Mki67", "Mpzl2", "1500009L16Rik"),
            reduction = "umap.harmony")

DimPlot(all_m_thymocyte, reduction = "umap.PCA", split.by = "genotype", group.by = "PCA_clusters")
FeaturePlot(all_m_thymocyte, 
            features = c("Cd24a", "Cd69", "Sell", "Cd4", "Cd8b1", "Foxp3", "Tcrg-C2", "Mki67", "Mpzl2", "1500009L16Rik"),
            reduction = "umap.PCA")

DimPlot(all_m_thymocyte, reduction = "umap.PCA", split.by = "genotype", group.by = "PCA_clusters")+
  xlim(-8, 3) + ylim(0, 15)

FeaturePlot(all_m_thymocyte, 
            features = c("CD4", "CD8b(Ly-3)", "CD8a"), split.by = "genotype",
            reduction = "umap.PCA")

Idents(all_m_thymocyte) <- all_m_thymocyte$PCA_clusters

# clustering
DimPlot(all_m_thymocyte, reduction = "umap.PCA", group.by = "PCA_clusters", label = T, label.box = T, repel = T)+NoLegend()
DotPlot(all_m_thymocyte, features = # PMID: 37580604 marker genes
          c("Snca", "Alas2", "Gypa", "Hba-a2", "Rsad2", "Oasl2", #Erythrocyte
            "Iglc3", "Iglc2", "Gm43291", "Ms4a1", "Ctsh", #B cells
            "Spi1", "Tyrobp", "Mpeg1", #Myeloid
            "Tcrg-C2", "Cpa3", "Gm16602", "Tcrg-C4", "5830411N06Rik", "Hes1", #GD T cells
            "Mpzl2", "1500009L16Rik", #DN
            "Hist1h2af", "Hist1h3b", "Esco2", "Mki67", "BC030867",#DP (P), Mature cycling
            "Insl5", "Slc7a11",#DP (Q1)
            "P2rx1", "Rag1", "BB031773", "Tctex1d1", #DP (Q2)
            "Rapsn", "Ccr4", #DP (Sig)
            "Ifit1", "Iigp1", #Interferon sig.
            "Itm2a", "Zbtb7b", "Cd40lg", "Chl1", "Gpr83", "Acvrl1",#Immature CD4, Neg. sel(2)
            "Gm26771", "Tesc", #Immature CD8
            "S1pr1", "Slfn1", #Mature CD4
            "Cd7",#NKT
            "Fam101b", "Dapl1", "Gm43698",#Mature CD8
            "Nkg7", #NKT
            "Egr2", "Nr4a3", "Pdcd1lg2", "2610204G07Rik", #Neg. sel. (1)
            "Tnfrsf4", "Tnfrsf9", "Arhgap20", "Foxp3", "Slc22a2", "Lad1", #Treg
            "S100a4", "Zbtb16", "Klrb1a", "Klra1" #NKT
          ), group.by =  "PCA_clusters")+RotatedAxis()+
  scale_size_area(
    max_size = 5,                         # visual max dot (adjust)
    limits   = c(0, 100),                 # clamp range
    breaks   = c(1, 10, 100), labels = c("1", "10", "100"), trans="log1p")

DimPlot(subset(all_m_thymocyte, subset= PCA_clusters %in% c("15", "12", "11", "18", "23")),reduction = "umap.PCA")

FeaturePlot(subset(all_m_thymocyte, subset= PCA_clusters %in% c("1", "3", "30", "31", "20", "23", "18", "7", "10", "8",
                                            "5", "19", "27", "11", "15", "12", "2", "9")),
            features = c(
              "Ikzf2", "Bcl2l11", # PMID: 23337809 genes
              "Nr4a1", "Nr4a2", "Nr4a3", # PMID: 39824797 genes
              "Pdcd1", "Il2rb", "Nrgn", "Cd160", "Trpm1", # PMID: 33242395 genes, IELp, Sig-4
              "Trib1", "Cd69", #  PMID: 33242395 genes, TCR induced, Sig-5a
              "Nfkb1", "Nfkb2", "Nfkbia", "Nfkbid", "Nfkbiz", # PMID: 33242395 genes, NFkB family, Sig-5a
              "Tnfrsf18", "Tnfrsf4" # PMID: 33242395 genes, pre-Treg, Sig-5b
            ),
            reduction = "umap.PCA")

DotPlot(subset(all_m_thymocyte, subset= PCA_clusters %in% c("1", "3", "30", "31", "20", "23", "18", "7", "10", "8",
                                        "5", "19", "27", "11", "15", "12", "2", "9")),
        features = c(
          "Ikzf2", "Bcl2l11", # PMID: 23337809 genes
          "Nr4a1", "Nr4a2", "Nr4a3", # PMID: 39824797 genes
          "Pdcd1", "Il2rb", "Nrgn", "Cd160", "Trpm1", # PMID: 33242395 genes, IELp, Sig-4
          "Trib1", "Cd69", #  PMID: 33242395 genes, TCR induced, Sig-5a
          "Nfkb1", "Nfkb2", "Nfkbia", "Nfkbid", "Nfkbiz", # PMID: 33242395 genes, NFkB family, Sig-5a
          "Tnfrsf18", "Tnfrsf4" # PMID: 33242395 genes, pre-Treg, Sig-5b
        ))+RotatedAxis()

DotPlot(all_m_thymocyte, features = c( # PMID: 23337809 genes
  "Il2ra", "Foxp3", "Cd24a", "Cd69", "Ikzf2", "Pdcd1",  "Bcl2l11", "Card11", "Rel", "Ccr7", "Cd4", "Cd8b1"),
  group.by =  "PCA_clusters")+RotatedAxis()

all_m_thymocyte <- RenameIdents(all_m_thymocyte, "0"= "DP", "1"="Pre DP-sig", "2"="Mature CD4-1", "3"="DP-sig" , "4"="DP" , "5"="Immature CD8-2", "6"="DP",
                                "7"="Immature CD8-1", "8"="Immature CD4-2", "9"="Mature CD4-2", "10"="Immature CD4-1", "11"="Mature CD8-1", "12"="Treg",
                                "13"="Erythrocyte", "14"="DP-P", "15"="intCD4", "16"="DN", "17"="gdT", "18"="AgoSel1",
                                "19"="Mature CD8-2", "20"="NKT", "21"="gdT", "22"="DP", "23"="AgoSel2", "24"="Erythrocyte",
                                "25"="DP-P", "26"="DP-P", "27"="Mature cycling T", "28"= "Myeloid", "29"="Myeloid", "30"="Pre DP-sig",
                                "31"="Pre DP-sig", "32"= "B cell")

DimPlot(all_m_thymocyte, reduction = "umap.PCA", label = TRUE, label.box = T, repel = T)

## add Pre-Treg cluster
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
          yaxis=list(title="SP_2")
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

treg_cells <- WhichCells(all_m_thymocyte, idents = "Treg")

coords <- Embeddings(all_m_thymocyte[["umap.PCA"]])[treg_cells, 1:2]

df <- data.frame(cell = treg_cells, x    = coords[,1], y    = coords[,2])

picked_cells <- lasso_select(df)

new_idents <- Idents(all_m_thymocyte)
levels(new_idents) <- c(levels(new_idents), "Pre-Treg")
new_idents[picked_cells] <- "Pre-Treg"
Idents(all_m_thymocyte) <- new_idents

DimPlot(all_m_thymocyte, reduction = "umap.PCA", label = TRUE, label.box = T, repel = T)

all_m_thymocyte <- SetIdent(all_m_thymocyte, value = factor(Idents(all_m_thymocyte), levels = c(
  "Erythrocyte","B cell", "Myeloid", "gdT",
  "DN", "DP-P", "DP", "Pre DP-sig", "DP-sig", 
  "Immature CD4-1",  "Immature CD4-2", "intCD4", "Mature CD4-1", "Mature CD4-2",  
  "Immature CD8-1", "Immature CD8-2", "Mature CD8-1", "Mature CD8-2", "Mature cycling T", 
  "AgoSel1", "Pre-Treg", "Treg", 
  "AgoSel2", "NKT"
)))

DimPlot(all_m_thymocyte, reduction = "umap.PCA", split.by = "genotype")

saveRDS(all_m_thymocyte, "PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")
all_m_thymocyte <- readRDS("PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")

#plotting
DimPlot(all_m_thymocyte, reduction = "umap.PCA", label = T, label.box = T, repel = T)
ggsave(filename = "pic/FigS1B_1.pdf", plot = get_last_plot(), width = 11, height = 7)

all_m_thymocyte$genotype_highlight <- ""
all_m_thymocyte$genotype_highlight[all_m_thymocyte$genotype 
                                   %in% c("MHCI_KO", "OT_II", "AND")] <- "CD4-fated"
all_m_thymocyte$genotype_highlight[all_m_thymocyte$genotype 
                                   %in% c("MHCII_KO", "OT_I", "F5")] <- "CD8-fated"
all_m_thymocyte$genotype_highlight[all_m_thymocyte$genotype 
                                   %in% c("WT")] <- "WT"

all_m_thymocyte$genotype_highlight <- factor(all_m_thymocyte$genotype_highlight,
                                             levels = c("CD4-fated", "CD8-fated", "WT"))

DimPlot(all_m_thymocyte, reduction = "umap.PCA", group.by  = "genotype_highlight",
        cols = c("CD4-fated"="green", "CD8-fated"="red","WT"="grey"))+
  ggtitle(NULL) & theme(axis.title = element_blank())
ggsave(filename = "pic/FigS1C.pdf", plot = get_last_plot(), width = 5, height = 4)


FeaturePlot(all_m_thymocyte, reduction = "umap.PCA", features = c("Zbtb7b", "Runx3"))&
  theme(axis.title = element_blank())
ggsave(filename = "pic/FigS1C_2.pdf", plot = get_last_plot(), width = 7.5, height = 4)

DotPlot(all_m_thymocyte, scale.by = "size", dot.scale = 6,
        features = c(
          "Hba-a1", "Hba-a2", "Gypa", #Erythrocyte
          "Iglc2", "Iglc3",  #B cells
          "Tyrobp", "Mpeg1", #Myeloid
          "Tcrg-C2", "Tcrg-C4", #GD T cells
          "Mpzl2", "1500009L16Rik", #DN
          "Mki67", "BC030867",#DP (P), Mature cycling
          "Slc7a11",#DP (Q1)
          "Rag1", #DP (Q2)
          "Rapsn", "Ccr4", #DP (Sig)
          "Itm2a", "Chl1", "Cd40lg", "Zbtb7b", #Immature CD4, Neg. sel(2)
          "S1pr1", "Slfn1", "Sell", #Mature CD4
          "Runx3","Gm26771",#Immature CD8
          "Fam101b", "Dapl1", #Mature CD8
          "Egr2", "Pdcd1lg2", "2610204G07Rik", #Neg. sel. (1)
          "Tnfrsf9", "Arhgap20", "Foxp3", #Treg
          "S100a4", "Zbtb16", #NKT
          "Ikzf2", "Bcl2l11", # PMID: 23337809 genes
          "Nr4a1", "Nr4a3",  # PMID: 39824797 genes
          "Trib1", "Cd69", #  PMID: 33242395 genes, TCR induced, Sig-5a
          "Nfkb1", "Nfkbia", # PMID: 33242395 genes, NFkB family, Sig-5a
          "Tnfrsf4", "Tnfrsf18",  # PMID: 33242395 genes, pre-Treg, Sig-5b
          "Pdcd1", "Nrgn", "Cd160" # PMID: 33242395 genes, IELp, Sig-4
          ))+RotatedAxis()+scale_y_discrete(limits = rev(levels(all_m_thymocyte)))+
  labs(x=NULL, y=NULL)
ggsave(filename = "pic/FigS1D.pdf", plot = get_last_plot(), width = 15, height = 6)

# Rat vs Mouse thymocyte integration
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

RatThymocyteSimple.obj <- subset(RatThymocyte.obj, ident=c("DN", "DP-P1", "DP-P2", "DP-P3", "DP1", "DP2", "Early SP", "Mature SP"))
Idents(RatThymocyteSimple.obj) <- factor(paste0("Rat: ", Idents(RatThymocyteSimple.obj)))

MouseThymocyteSimple.obj <-  subset(all_m_thymocyte, ident= c("DN", "DP-P", "DP", "Pre DP-sig", "DP-sig", 
                                                              "Immature CD4-1", "Immature CD4-2", "intCD4", "Mature CD4-1", "Mature CD4-2", 
                                                              "Immature CD8-1", "Immature CD8-2", "Mature CD8-1", "Mature CD8-2", "Mature cycling T",
                                                              "AgoSel1", "Pre-Treg", "Treg", "AgoSel2", "NKT", "gdT"))

Idents(MouseThymocyteSimple.obj) <- factor(paste0("Mouse: ", Idents(MouseThymocyteSimple.obj)))

MouseThymocyteSimple.obj$orig.ident <- "MouseThymocyte"

class(MouseThymocyteSimple.obj[["RNA"]])
MouseThymocyteSimple.obj <- UpdateSeuratObject(MouseThymocyteSimple.obj)
MouseThymocyteSimple.obj[["RNA"]] <- as(object = MouseThymocyteSimple.obj[["RNA"]], Class = "Assay5")
class(MouseThymocyteSimple.obj[["RNA"]])

library(orthogene)
rownames(MouseThymocyteSimple.obj) <- make.unique(convert_orthologs(rownames(MouseThymocyteSimple.obj), gene_input = "rownames", 
                                                                 gene_output = "columns", input_species = "mouse",
                                                                 output_species = "rat", drop_nonorths = F,
                                                                 non121_strategy = 5)[,2], sep = "_")

Rat_MouseThymocyte.obj <- merge(RatThymocyteSimple.obj, y= MouseThymocyteSimple.obj, add.cell.ids = c("RatThymocyte", "MouseThymocyte"))

Rat_MouseThymocyte.obj$merged_cluster <- Idents(Rat_MouseThymocyte.obj)

Rat_MouseThymocyte.obj <- Rat_MouseThymocyte.obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
Rat_MouseThymocyte.obj <- NormalizeData(Rat_MouseThymocyte.obj, normalization.method = "CLR", margin = 2, assay = "ADT")

DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.unintegrated", split.by = "orig.ident")

Rat_MouseThymocyte.obj <- IntegrateLayers(Rat_MouseThymocyte.obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
Rat_MouseThymocyte.obj <- FindNeighbors(Rat_MouseThymocyte.obj, reduction = "integrated.cca", dims = 1:30)
Rat_MouseThymocyte.obj <- FindClusters(Rat_MouseThymocyte.obj, resolution = 0.4, cluster.name = "cca_clusters")
Rat_MouseThymocyte.obj  <- RunUMAP(Rat_MouseThymocyte.obj , reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.cca", split.by = "orig.ident", label = T, label.box = T, repel = T)

Rat_MouseThymocyte.obj <- IntegrateLayers(Rat_MouseThymocyte.obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", assay = "RNA")
Rat_MouseThymocyte.obj <- FindNeighbors(Rat_MouseThymocyte.obj, reduction = "integrated.rpca", dims = 1:30,assay = "RNA")
Rat_MouseThymocyte.obj <- FindClusters(Rat_MouseThymocyte.obj, resolution = 0.4, cluster.name = "rpca_clusters", assay = "RNA")
Rat_MouseThymocyte.obj  <- RunUMAP(Rat_MouseThymocyte.obj , reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", assay = "RNA")
DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.rpca", split.by ="orig.ident", label = T, label.box = T)

Rat_MouseThymocyte.obj <- IntegrateLayers(Rat_MouseThymocyte.obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
Rat_MouseThymocyte.obj <- FindNeighbors(Rat_MouseThymocyte.obj, reduction = "harmony", dims = 1:30)
Rat_MouseThymocyte.obj <- FindClusters(Rat_MouseThymocyte.obj, resolution = 0.4, cluster.name = "harmony_clusters")
Rat_MouseThymocyte.obj  <- RunUMAP(Rat_MouseThymocyte.obj , reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.harmony", split.by ="orig.ident", label = T, label.box = T)

Rat_MouseThymocyte.obj <- IntegrateLayers(Rat_MouseThymocyte.obj, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.PCA")
Rat_MouseThymocyte.obj <- FindNeighbors(Rat_MouseThymocyte.obj, reduction = "integrated.PCA", dims = 1:30)
Rat_MouseThymocyte.obj <- FindClusters(Rat_MouseThymocyte.obj, resolution = 0.4, cluster.name = "PCA_clusters")
Rat_MouseThymocyte.obj  <- RunUMAP(Rat_MouseThymocyte.obj , reduction = "integrated.PCA", dims = 1:30, reduction.name = "umap.PCA")
DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.PCA", split.by ="orig.ident", label = T, label.box = T)

Rat_MouseThymocyte.obj <- SetIdent(Rat_MouseThymocyte.obj, value = Rat_MouseThymocyte.obj$merged_cluster)

Rat_MouseThymocyte.obj<-SetIdent(Rat_MouseThymocyte.obj, 
                                 value =factor(Idents(Rat_MouseThymocyte.obj), 
                                               levels = c("Mouse: gdT", "Mouse: DN", "Mouse: DP-P", "Mouse: DP", 
                                                          "Mouse: Pre DP-sig", "Mouse: DP-sig", 
                                                          "Mouse: Immature CD4-1", "Mouse: Immature CD4-2", "Mouse: intCD4", 
                                                          "Mouse: Mature CD4-1", "Mouse: Mature CD4-2", 
                                                          "Mouse: Immature CD8-1", "Mouse: Immature CD8-2", 
                                                          "Mouse: Mature CD8-1", "Mouse: Mature CD8-2", "Mouse: Mature cycling T", 
                                                          "Mouse: AgoSel1", "Mouse: Pre-Treg", "Mouse: Treg",
                                                          "Mouse: AgoSel2", "Mouse: NKT", 
                                                          "Rat: DN", "Rat: DP-P1", "Rat: DP-P2", "Rat: DP-P3",
                                                          "Rat: DP1", "Rat: DP2", "Rat: Early SP", "Rat: Mature SP")))

saveRDS(Rat_MouseThymocyte.obj, "20250904 Rat thymocyte scRNA-seq/Rat_MouseThymocyte.obj.rds")
Rat_MouseThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/Rat_MouseThymocyte.obj.rds")

Rat_MouseThymocyte.obj$orig.ident <- factor(Rat_MouseThymocyte.obj$orig.ident, levels = c("RatThymocyte", "MouseThymocyte"))
DimPlot(Rat_MouseThymocyte.obj, reduction = "umap.PCA", split.by = "orig.ident", 
        label = T, label.box = T, repel = T, label.size = 3)+NoLegend()
ggsave(filename = "pic/FigS1E.pdf", plot = get_last_plot(), width = 10, height = 5)

# Project mouse subsets to rat dataset
project_labels_from_mouse <- function(
    obj,
    reduction     = c("cca","pca","umap.cca","umap"),
    dims          = 1:30,
    k             = 50,
    vote_min_frac = 0.55,             # majority needed to accept the label
    set_idents    = TRUE,             # write back to Idents()
    projected_prefix = "Projected: ", # final label prefix
    mouse_label_prefix_regex = "^Mouse:\\s*" # remove this from donor labels
) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Please install 'FNN'.")
  
  # pick first available reduction from choices
  reduction <- reduction[reduction %in% names(obj@reductions)][1]
  if (is.na(reduction)) stop("Chosen reduction not found in object.")
  
  emb_all <- Seurat::Embeddings(obj, reduction = reduction)
  
  # make sure dims exist
  dims <- dims[dims <= ncol(emb_all)]
  if (!length(dims)) stop("No valid dimensions found in the chosen reduction.")
  emb_all <- emb_all[, dims, drop = FALSE]
  
  # identify species by orig.ident
  stopifnot("orig.ident" %in% colnames(obj@meta.data))
  is_mouse <- obj$orig.ident == "MouseThymocyte"
  is_rat   <- obj$orig.ident == "RatThymocyte"
  if (!any(is_mouse) || !any(is_rat)) stop("No mouse or rat cells found via orig.ident.")
  
  mouse_cells <- which(is_mouse)
  rat_cells   <- which(is_rat)
  
  # donor labels (fine-grained) from current Idents
  mouse_labels <- as.character(Seurat::Idents(obj))[mouse_cells]
  
  # KNN from each rat cell to mouse cells
  nn <- FNN::get.knnx(
    data  = emb_all[mouse_cells, , drop = FALSE],
    query = emb_all[rat_cells,   , drop = FALSE],
    k     = k
  )$nn.index
  
  # majority vote per rat cell -> "Projected: XX"
  rat_proj <- vapply(seq_len(nrow(nn)), function(i) {
    labs <- mouse_labels[nn[i, ]]
    tab  <- sort(table(labs), decreasing = TRUE)
    top  <- names(tab)[1]
    frac <- tab[1] / sum(tab)
    if (frac >= vote_min_frac) {
      # strip any "Mouse: " (or your regex) then add "Projected: "
      top_clean <- sub(mouse_label_prefix_regex, "", top)
      paste0(projected_prefix, top_clean)
    } else NA_character_
  }, character(1))
  
  # store to meta
  out <- rep(NA_character_, ncol(obj)); names(out) <- colnames(obj)
  out[rat_cells] <- rat_proj
  obj$projected_label <- out
  
  # optionally set Idents: keep mouse as-is; replace rat where we have a projection
  if (set_idents) {
    cur <- as.character(Seurat::Idents(obj))
    cur[rat_cells] <- ifelse(is.na(rat_proj), cur[rat_cells], rat_proj)
    Seurat::Idents(obj) <- factor(cur)
  }
  
  obj
}

Rat_MouseThymocyte.obj <- project_labels_from_mouse(
  Rat_MouseThymocyte.obj,
  reduction     = "umap.PCA",
  dims          = 1:50,
  k             = 5,
  vote_min_frac = 0.2,
  set_idents    = TRUE
)

table(Rat_MouseThymocyte.obj$projected_label)

Rat_MouseThymocyte.obj$projected_label <- factor(Rat_MouseThymocyte.obj$projected_label,
                                                 levels = c("Projected: DN", "Projected: DP-P", "Projected: DP",
                                                            "Projected: Pre DP-sig",  "Projected: DP-sig", 
                                                            "Projected: Immature CD4-1", "Projected: Immature CD4-2", "Projected: intCD4",
                                                            "Projected: Mature CD4-1", "Projected: Mature CD4-2",   
                                                            "Projected: Immature CD8-1", "Projected: Immature CD8-2", 
                                                            "Projected: Mature CD8-1", "Projected: Mature CD8-2", "Projected: Mature cycling T", 
                                                            "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg", 
                                                            "Projected: AgoSel2", "Projected: NKT", "Projected: gdT"))

DimPlot(subset(Rat_MouseThymocyte.obj, subset=orig.ident=="RatThymocyte"), 
        reduction = "umap.PCA", group.by = "projected_label", 
        label = T, label.box = T, repel = T, label.size = 3)+NoLegend() +
  ggtitle("Projected labels to Rat")
ggsave(filename = "pic/FigS1E_2.pdf", plot = get_last_plot(), width = 5, height = 5)

# Add projected clustering to the original rat thymocyte object
proj <- Rat_MouseThymocyte.obj$projected_label
names(proj) <- colnames(Rat_MouseThymocyte.obj)
proj <- proj[grep("^RatThymocyte_", names(proj))]
names(proj) <- sub("^RatThymocyte_", "", names(proj))
aligned <- proj[colnames(RatThymocyte.obj)]
names(aligned) <- colnames(RatThymocyte.obj)
RatThymocyte.obj <- Seurat::AddMetaData(RatThymocyte.obj, metadata = aligned, col.name = "projected_label")

RatThymocyte.obj$projected_label <- factor(RatThymocyte.obj$projected_label, 
                                           levels = c("Projected: DN", "Projected: DP-P", "Projected: DP",
                                                      "Projected: Pre DP-sig",  "Projected: DP-sig", 
                                                      "Projected: Immature CD4-1", "Projected: Immature CD4-2", "Projected: intCD4",
                                                      "Projected: Mature CD4-1", "Projected: Mature CD4-2",   
                                                      "Projected: Immature CD8-1", "Projected: Immature CD8-2", 
                                                      "Projected: Mature CD8-1", "Projected: Mature CD8-2", "Projected: Mature cycling T", 
                                                      "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg", 
                                                      "Projected: AgoSel2", "Projected: NKT", "Projected: gdT"
                                                      ))

DimPlot(RatThymocyte.obj, group.by = "projected_label")

saveRDS(RatThymocyte.obj, "20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

# GSEA between rat mature and immature thymocytes
RatThy_MatureImmature_markers <- FindMarkers(Seurat::SetIdent(RatThymocyte.obj, value=RatThymocyte.obj$projected_label),
                                            ident.1 = c("Projected: Mature CD4-1", "Projected: Mature CD4-2",
                                                        "Projected: Mature CD8-1", "Projected: Mature CD8-2",
                                                        "Projected: Pre Treg", "Projected: Treg"),
                                            ident.2 = c("Projected: DP-sig", "Projected: Immature CD4-1", 
                                                        "Projected: Immature CD8-1", "Projected: AgoSel1"),
                                            only.pos=F, min.pct=0.1, logfc.threshold = 0.1)

RatThy_MatureImmature_markers <- RatThy_MatureImmature_markers[order(RatThy_MatureImmature_markers$avg_log2FC, decreasing=T),]
RatThy_MatureImmature_markers$genes <- rownames(RatThy_MatureImmature_markers)
View(RatThy_MatureImmature_markers)

RatThy_MatureImmature_list <- RatThy_MatureImmature_markers$avg_log2FC
names(RatThy_MatureImmature_list) <- rownames(RatThy_MatureImmature_markers)

RatThy_MatureImmature_list <- sort(RatThy_MatureImmature_list, decreasing = T)

RatThy_MatureImmature_GESA <- gseGO(geneList =RatThy_MatureImmature_list, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL", pAdjustMethod="none")

View(RatThy_MatureImmature_GESA@result)

RatThy_Mature_GESA <- RatThy_MatureImmature_GESA; RatThy_Immature_GESA <- RatThy_MatureImmature_GESA
RatThy_Mature_GESA@result <- RatThy_Mature_GESA@result[RatThy_Mature_GESA@result$enrichmentScore > 0,]
RatThy_Immature_GESA@result <- RatThy_Immature_GESA@result[RatThy_Immature_GESA@result$enrichmentScore < 0,]

treeplot(pairwise_termsim(RatThy_Mature_GESA)) | treeplot(pairwise_termsim(RatThy_Immature_GESA))
gseaplot2(RatThy_MatureImmature_GESA, geneSetID = 1:3)

# Rat thymic stroma object----
RatThymicStroma.obj <- CreateSeuratObject(counts = Read10X(
  "20240331 rat thymic stroma scPCR recounted/filtered_feature_bc_matrix/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
), project = "RatThymicStroma")
RatThymicStroma.obj

RatThymicStroma.obj[["percent.mt"]] <- PercentageFeatureSet(RatThymicStroma.obj, pattern = "^mt-")
head(RatThymicStroma.obj@meta.data, 5)
VlnPlot(RatThymicStroma.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mean(RatThymicStroma.obj@meta.data$nFeature_RNA)

plot1 <- FeatureScatter(RatThymicStroma.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RatThymicStroma.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

RatThymicStroma.obj <- subset(RatThymicStroma.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 5)

RatThymicStroma.obj <- NormalizeData(RatThymicStroma.obj)

RatThymicStroma.obj <- FindVariableFeatures(RatThymicStroma.obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(RatThymicStroma.obj), 10)

plot1 <- VariableFeaturePlot(RatThymicStroma.obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(RatThymicStroma.obj)
RatThymicStroma.obj <- ScaleData(RatThymicStroma.obj, features = all.genes)


RatThymicStroma.obj <- RunPCA(RatThymicStroma.obj, features = VariableFeatures(object = RatThymicStroma.obj))
print(RatThymicStroma.obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(RatThymicStroma.obj, dims = 1:2, reduction = "pca")
DimPlot(RatThymicStroma.obj, reduction = "pca")

DimHeatmap(RatThymicStroma.obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(RatThymicStroma.obj, dims = 1:15, cells = 500, balanced = TRUE)

RatThymicStroma.obj <- JackStraw(RatThymicStroma.obj, num.replicate = 100)
RatThymicStroma.obj <- ScoreJackStraw(RatThymicStroma.obj, dims = 1:15)
JackStrawPlot(RatThymicStroma.obj, dims = 1:15)
ElbowPlot(RatThymicStroma.obj)

RatThymicStroma.obj <- FindNeighbors(RatThymicStroma.obj)
RatThymicStroma.obj <- FindClusters(RatThymicStroma.obj, resolution = 0.5)
head(Idents(RatThymicStroma.obj), 5)

RatThymicStroma.obj <- RunUMAP (RatThymicStroma.obj, dims = 1:17)
DimPlot(RatThymicStroma.obj, label = TRUE, reduction = "umap")

cluster_names <- levels(Idents(RatThymicStroma.obj))

# Cluster identification
DotPlot(
  RatThymicStroma.obj,
  features = c("Ptprc", "Cd3d", "Trb", "Cd19", "Cd79a", "Ighm", "Itgam", "Klrb1a", "S100a8", "Foxn1","Epcam",
               "Krt8", "Ly75", "Aire", "Dcn", "Col3a1", "Twist2", "Pdpn", "Pdgfra", "Pdgfrb", "Dpp4" , "Pi16",
               "Mfap5", "Pecam1", "Ramp2", "Cav1", "Plvap", "Msln", "Lrrn4", "Upk3b", "Des", "Cox4i2", "Tagln", "Acta2", "Cspg4", "Itga7")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#name the clusters
ThymicStromaClusterNames <- c("TMC2", "TEC-b", "TMC1", "capFb", "TMC4", "TMC3", "Endo-1", "Myeloid", "vSMC/PC", 
                              "Thymocyte", "TEC-d", "TEC-c", "Endo-2", "TEC-a", "B cell", "TEC-e")
names(ThymicStromaClusterNames) <- levels(Idents(RatThymicStroma.obj))
RatThymicStroma.obj <- RenameIdents(object = RatThymicStroma.obj, ThymicStromaClusterNames)

#set order of clusters
ThymicCluster_order <- c("Thymocyte", "B cell", "Myeloid", 
                         "capFb", "TMC1", "TMC2", "TMC3", "TMC4", "vSMC/PC",
                         "Endo-1", "Endo-2", 
                         "TEC-a", "TEC-b", "TEC-c", "TEC-d", "TEC-e")
RatThymicStroma.obj <- SetIdent(RatThymicStroma.obj, value = factor(Idents(RatThymicStroma.obj), levels = ThymicCluster_order))

DimPlot(RatThymicStroma.obj, label = TRUE, label.box = T, repel = T, reduction = "umap")
ggsave(filename = "pic/FigS2A.pdf", plot = get_last_plot(), width = 7, height = 5)

DotPlot(
  RatThymicStroma.obj,
  features = c("Ptprc", "Cd3d", "Cd19", "Cd79a", "Itgam", "S100a8", "Dpp4" , "Pi16", "Mfap5",
               "Dcn", "Col3a1", "Twist2", "Pdpn", "Pdgfra", "Pdgfrb",
               "Des", "Cox4i2", "Tagln", "Acta2","Itga7",
               "Pecam1", "Ramp2", "Cav1", "Plvap", 
               "Foxn1", "Krt5", "Krt8", "Epcam","Ly75", "Aire")) + 
  scale_y_discrete(limits = rev(levels(RatThymicStroma.obj))) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(filename = "pic/FigS2A_2.pdf", plot = get_last_plot(), width = 10, height = 4.2)

# save/load
saveRDS(RatThymicStroma.obj, file = "20240331 rat thymic stroma scPCR recounted/RatThymicStroma.rds")
RatThymicStroma.obj <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicStroma.rds")

# Read mouse stromal cells and project labels to rat stromal cells----
library(zellkonverter)

sce <- readH5AD(
  "PMID_39112630 mouse human integrated thymic stroma/thymosight_cd45neg_TOTAL_mouse.h5ad",
  use_hdf5 = F
)

cd <- as.data.frame(colData(sce))
cd[] <- lapply(cd, function(x) if (is.factor(x)) as.character(x) else x)
colData(sce) <- S4Vectors::DataFrame(cd, row.names = colnames(sce))

MouseThymicStroma.obj <- as.Seurat(sce, counts = "X", data =  NULL)

MouseThymicStroma.obj <- RenameAssays(MouseThymicStroma.obj, "originalexp" = "RNA")
DefaultAssay(MouseThymicStroma) <- "RNA"

Idents(MouseThymicStroma.obj) <- "cell_type_subset"

DimPlot(MouseThymicStroma.obj, label = T, label.box = T, repel = T)+NoLegend()

# Save
saveRDS(MouseThymicStroma.obj, "PMID_39112630 mouse human integrated thymic stroma/MouseThymicStroma.rds")
MouseThymicStroma.obj <- readRDS(file = "PMID_39112630 mouse human integrated thymic stroma/MouseThymicStroma.rds")

#downsample the object
set.seed(42)
idents <- Idents(MouseThymicStroma.obj)
cell_counts <- table(idents)

base_target <- round(as.numeric(cell_counts) / sum(cell_counts) * 40000) # base target (proportional to 20,000)
names(base_target) <- names(cell_counts)

orig_n <- as.numeric(cell_counts)

target_per_group <- ifelse( # rule: - if a cluster has <200 cells total -> keep all, - else keep max(200, base_target)
  orig_n < 200,
  orig_n,
  pmax(200, base_target)
)

target_per_group <- pmin(target_per_group, orig_n) # never ask for more than exist
names(target_per_group) <- names(cell_counts)

sampled_cells <- unlist(lapply(names(target_per_group), function(g) { # sample
  cells_g <- WhichCells(MouseThymicStroma.obj, idents = g)
  k <- target_per_group[[g]]
  if (k >= length(cells_g)) cells_g else sample(cells_g, k)
}), use.names = FALSE)

MouseThymicStroma_l.obj <- subset(MouseThymicStroma.obj, cells = sampled_cells)

DimPlot(MouseThymicStroma_l.obj, label = TRUE, label.box = TRUE, repel = TRUE) + NoLegend()

#convert_orthologs
class(MouseThymicStroma_l.obj[["RNA"]])

MouseThymicStroma_l_ortho.obj <- MouseThymicStroma_l.obj
MouseThymicStroma_l_ortho.obj <- UpdateSeuratObject(MouseThymicStroma_l_ortho.obj)
MouseThymicStroma_l_ortho.obj[["RNA"]] <- as(object = MouseThymicStroma_l_ortho.obj[["RNA"]], Class = "Assay5")
rownames(MouseThymicStroma_l_ortho.obj) <- make.unique(convert_orthologs(rownames(MouseThymicStroma_l_ortho.obj), gene_input = "rownames", 
                                                                       gene_output = "columns", input_species = "mouse",
                                                                       output_species = "rat", drop_nonorths = F,
                                                                       non121_strategy = 5)[,2], sep = "_")
FeaturePlot(MouseThymicStroma_l.obj, "Ccl21a")

FeaturePlot(MouseThymicStroma_l_ortho.obj, "Ccl21")

#integrate mouse and rat stroma
RatThymicStroma.obj <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicStroma.rds")

RatThymicStroma_c.obj <- subset(RatThymicStroma.obj, idents= setdiff(levels(RatThymicStroma.obj), c("Myeloid", "B cell", "Thymocyte")))
Idents(RatThymicStroma_c.obj) <- paste0("Rat: ", Idents(RatThymicStroma_c.obj))


#MouseThymicStroma_l_ortho.obj <- subset(MouseThymicStroma_l_ortho.obj, 
#                                        idents= setdiff(levels(MouseThymicStroma_l_ortho.obj), c("FB:fetal", "TEC:aaTEC1", "TEC:aaTEC2", "vSMC/PC:fetal", "EC:fetal")))

Idents(MouseThymicStroma_l_ortho.obj) <- paste0("Mouse: ", Idents(MouseThymicStroma_l_ortho.obj))

RatMouseStroma.obj <- merge(RatThymicStroma_c.obj, y = MouseThymicStroma_l_ortho.obj,
                            add.cell.ids = c("RatStroma", "MouseStroma"), project="RatMouseThymicStroma")

RatMouseStroma.obj$merged_cluster <- Idents(RatMouseStroma.obj)

RatMouseStroma.obj <- RatMouseStroma.obj %>% NormalizeData()%>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  FindNeighbors(dims=1:20, reduction = "pca") %>% FindClusters(resolution = 0.3, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

RatMouseStroma.obj$species <- ifelse(grepl("^Rat", RatMouseStroma.obj@meta.data$merged_cluster), "Rat", "Mouse")

DimPlot(RatMouseStroma.obj, reduction = "umap.unintegrated", split.by = "species", group.by = "unintegrated_clusters", label = T)

DimPlot(RatMouseStroma.obj, reduction = "umap.unintegrated", split.by = "species", group.by = "ident", label = T, label.box = T, repel = T)+NoLegend()

DimPlot(RatMouseStroma.obj, reduction = "pca", group.by = "ident", label = T, label.box = T, repel = T)+NoLegend()

RatMouseStroma.obj <- IntegrateLayers(object = RatMouseStroma.obj, method = CCAIntegration, orig.reduction = "pca",
                                  new.reduction ="integrated.cca")
RatMouseStroma.obj <- RatMouseStroma.obj %>% FindNeighbors(reduction = "integrated.cca", dims = 1:50) %>%
  FindClusters(resolution = 0.3,cluster.name = "cca_clusters") %>% 
  RunUMAP(reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")
DimPlot(RatMouseStroma.obj, reduction = "umap.cca",  label = T, label.box = T, repel = T, split.by = "species")+NoLegend()

RatMouseStroma.obj <- IntegrateLayers(object = RatMouseStroma.obj, method = RPCAIntegration, orig.reduction = "pca",
                                      new.reduction ="integrated.rpca")
RatMouseStroma.obj <- RatMouseStroma.obj %>% FindNeighbors(reduction = "integrated.rpca", dims = 1:50) %>%
  FindClusters(resolution = 0.3,cluster.name = "rpca_clusters") %>% 
  RunUMAP(reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca")
DimPlot(RatMouseStroma.obj, reduction = "umap.rpca",  label = T, label.box = T, repel = T, split.by = "species")+NoLegend()

RatMouseStroma.obj <- IntegrateLayers(object = RatMouseStroma.obj, method = HarmonyIntegration, orig.reduction = "pca",
                                      new.reduction ="integrated.harmony")
RatMouseStroma.obj <- RatMouseStroma.obj %>% FindNeighbors(reduction = "integrated.harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.3,cluster.name = "rpca_clusters") %>% 
  RunUMAP(reduction = "integrated.harmony", dims = 1:50, reduction.name = "umap.harmony")
DimPlot(RatMouseStroma.obj, reduction = "umap.harmony",  label = T, label.box = T, repel = T, split.by = "species")+NoLegend()

RatMouseStroma.obj <- IntegrateLayers(object = RatMouseStroma.obj, method = JointPCAIntegration, orig.reduction = "pca",
                                      new.reduction ="integrated.PCA")
RatMouseStroma.obj <- RatMouseStroma.obj %>% FindNeighbors(reduction = "integrated.PCA", dims = 1:50) %>%
  FindClusters(resolution = 0.3,cluster.name = "PCA_clusters") %>% 
  RunUMAP(reduction = "integrated.PCA", dims = 1:50, reduction.name = "umap.PCA")
DimPlot(RatMouseStroma.obj, reduction = "umap.PCA",  label = T, label.box = T, repel = T, split.by = "species")+NoLegend()

RatMouseStroma.obj <- SetIdent(RatMouseStroma.obj, value = factor(RatMouseStroma.obj$merged_cluster,
                               levels=c("Rat: capFb", "Rat: TMC1", "Rat: TMC2", "Rat: TMC3", "Rat: TMC4", "Rat: vSMC/PC",  
                                        "Rat: Endo-1", "Rat: Endo-2", 
                                        "Rat: TEC-a", "Rat: TEC-b", "Rat: TEC-c", "Rat: TEC-d", "Rat: TEC-e", 
                                        "Mouse: FB:capsFB", "Mouse: FB:intFB", "Mouse: FB:medFB",  "Mouse: vSMC/PC", "Mouse: vSMC/PC:fetal", "Mouse: Fat", "Mouse: FB:fetal", 
                                        "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC", "Mouse: MEC", "Mouse: EC:fetal",
                                        "Mouse: TEPC", "Mouse: TEC:early Pr", "Mouse: TEC:cTEC", "Mouse: TEC:mTEC1", "Mouse: TEC:mTEC-prol", "Mouse: TEC:mTEC2",
                                        "Mouse: TEC:aaTEC1", "Mouse: TEC:aaTEC2", 
                                        "Mouse: TEC:mimetic(tuft)", "Mouse: TEC:mimetic(parathyroid)", "Mouse: TEC:mimetic(basal)", 
                                        "Mouse: TEC:mimetic(microfold)", "Mouse: TEC:mimetic(goblet)", "Mouse: TEC:mimetic(muscle)", "Mouse: TEC:mimetic(neuroendo)",
                                        "Mouse: TEC:mimetic(ciliated)", "Mouse: TEC:mimetic(ionocyte)", "Mouse: nmSC"
                                        )))

saveRDS(RatMouseStroma.obj, "20240331 rat thymic stroma scPCR recounted/RatMouseStroma.rds")
RatMouseStroma.obj <- readRDS("20240331 rat thymic stroma scPCR recounted/RatMouseStroma.rds")

DimPlot(RatMouseStroma.obj, reduction = "umap.cca", group.by = "species") # check the distribution

RatMouseStroma.obj$species <- factor(RatMouseStroma.obj$species, levels = c("Rat", "Mouse"))

DimPlot(RatMouseStroma.obj, reduction = "umap.cca", label.size = 2,
        label = T, label.box = T, repel = T, split.by = "species")+NoLegend()
ggsave(filename = "pic/FigS2B.pdf", plot = get_last_plot(), width = 10, height = 5)

# Epithelial markers
DotPlot(RatMouseStroma.obj,
        features = c("Epcam", "Krt5", "Krt8", "Krt18", "Cldn3", "Cldn4", "Dsp",
                     "Aire", "Fezf2", "Tnfrsf11a", "Relb", "Ccl21"))+
  RotatedAxis() & theme(axis.title = element_blank())

FeaturePlot(RatMouseStroma.obj,
            reduction = "umap.cca", split.by = "species", 
            features = c("Aire", "Fezf2", "Tnfrsf11a", "Krt5", "Relb", "Ccl21")) & theme(axis.title = element_blank())

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj, 
                                     features = list(c("Epcam", "Krt5", "Krt8", "Krt18", "Cldn3", "Cldn4", "Dsp",
                                                       "Aire", "Fezf2", "Tnfrsf11a", "Krt5", "Relb", "Ccl21")),
                                     name = "Epithelial Score")
FeaturePlot(RatMouseStroma.obj,"Epithelial Score1", reduction = "umap.cca", split.by = "species")
VlnPlot(RatMouseStroma.obj,"Epithelial Score1")+NoLegend()

# GOBP Epithelial Differentiation score
library(msigdbr)
GOBP_df <- msigdbr(species = "Mus musculus", category = "C5")
r_ECD_genes <- GOBP_df[GOBP_df$gs_name=="GOBP_EPITHELIAL_CELL_DIFFERENTIATION",]$gene_symbol
r_ECD_genes <- convert_orthologs(unique(r_ECD_genes), input_species = "mouse", output_species = "rat", drop_nonorths = F,
                                 non121_strategy = 5, gene_output = "column")[, "ortholog_gene"]

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj, features = list(r_ECD_genes), name = "Epithelial cell differentiation")
FeaturePlot(RatMouseStroma.obj,"Epithelial cell differentiation1", reduction = "umap.cca", split.by = "species")
VlnPlot(RatMouseStroma.obj,"Epithelial cell differentiation1")+NoLegend()

# Mesenchymal markers
DotPlot(RatMouseStroma.obj,
        features = c("Vim", "Fn1", "Col1a1", "Col1a2", "Col4a1",
                     "Tagln", "Lox", "Sparc"))+
  RotatedAxis() & theme(axis.title = element_blank())

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj, 
                                     features = list(c("Vim", "Fn1", "Col1a1", "Col1a2", "Col4a1", "Tagln", "Lox", "Sparc",
                                                       "Postn", "Dcn", "Lum", "Pdgfra", "Pdgfrb")),
                                     name = "Mesenchymal Score")
FeaturePlot(RatMouseStroma.obj,"Mesenchymal Score1", reduction = "umap.cca", split.by = "species")
VlnPlot(RatMouseStroma.obj,"Mesenchymal Score1")+NoLegend()

# GOBP Mesenchymal cell Differentiation score
r_MCD_genes <- GOBP_df[GOBP_df$gs_name=="GOBP_MESENCHYMAL_CELL_DIFFERENTIATION",]$gene_symbol
r_MCD_genes <- convert_orthologs(unique(r_MCD_genes), input_species = "mouse", output_species = "rat", drop_nonorths = F,
                                 non121_strategy = 5, gene_output = "column")[, "ortholog_gene"]

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj, features = list(r_MCD_genes), name = "Mesenchymal cell differentiation")
FeaturePlot(RatMouseStroma.obj,"Mesenchymal cell differentiation1", reduction = "umap.cca", split.by = "species")
VlnPlot(RatMouseStroma.obj,"Mesenchymal cell differentiation1")+NoLegend()

# cell cycle
s.genes <- tools::toTitleCase(tolower(cc.genes$s.genes))
g2m.genes <- tools::toTitleCase(tolower(cc.genes$g2m.genes))

RatMouseStroma.obj <-CellCycleScoring(RatMouseStroma.obj, s.features = as.character(convert_orthologs(as.data.frame(s.genes), input_species = "mouse", output_species = "rat", non121_strategy = 5, gene_input = "s.genes", gene_output = "dict"))
                                                 , g2m.features = as.character(convert_orthologs(as.data.frame(g2m.genes), input_species = "mouse", output_species = "rat", non121_strategy = 5, gene_input = "g2m.genes", gene_output = "dict")))
VlnPlot(RatMouseStroma.obj, features = c("S.Score","G2M.Score"))

# EMT markers
FeaturePlot(RatMouseStroma.obj,
            reduction = "umap.cca", split.by = "species", 
            features = c("Zeb1", "Zeb2", "Snai1", "Twist1")) & theme(axis.title = element_blank())

# EMT gene score
library(msigdbr)

H_df <- msigdbr(species = "Mus musculus", category = "H")
H_df["r_gene_symbol"] <- convert_orthologs(H_df, input_species = "mouse", output_species = "rat", drop_nonorths = F,
                                           non121_strategy = 5, gene_input = "gene_symbol", gene_output = "column")[, "ortholog_gene"]
H_df <- rename(H_df, m_gene_symbol = gene_symbol)
H_df <- rename(H_df, gene_symbol = r_gene_symbol)

r_EMT_genes <- H_df[H_df$gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",]$gene_symbol

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj, features = list(r_EMT_genes), name = "EMT hallmark")
FeaturePlot(RatMouseStroma.obj,"EMT hallmark1", reduction = "umap.cca", split.by = "species")
VlnPlot(RatMouseStroma.obj,"EMT hallmark1")+NoLegend()

# Project mouse TEC/EC subsets to rat subsets
RatMouseTEC_EC.obj <- 
  subset(RatMouseStroma.obj, idents= c("Rat: Endo-1", 
                                       "Rat: Endo-2", "Rat: TEC-a", "Rat: TEC-b", "Rat: TEC-c", "Rat: TEC-d", 
                                       "Rat: TEC-e", 
                                       "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC", "Mouse: MEC", 
                                       "Mouse: EC:fetal", "Mouse: TEPC", "Mouse: TEC:early Pr", "Mouse: TEC:cTEC", 
                                       "Mouse: TEC:mTEC1", "Mouse: TEC:mTEC-prol", "Mouse: TEC:mTEC2", 
                                       "Mouse: TEC:aaTEC1", "Mouse: TEC:aaTEC2", "Mouse: TEC:mimetic(tuft)", 
                                       "Mouse: TEC:mimetic(parathyroid)", "Mouse: TEC:mimetic(basal)", 
                                       "Mouse: TEC:mimetic(microfold)", "Mouse: TEC:mimetic(goblet)", 
                                       "Mouse: TEC:mimetic(muscle)", "Mouse: TEC:mimetic(neuroendo)", 
                                       "Mouse: TEC:mimetic(ciliated)", "Mouse: TEC:mimetic(ionocyte)", 
                                       "Mouse: nmSC"))
  
project_labels_from_mouse <- function(
    obj,
    reduction     = c("cca","pca","umap.cca","umap"),
    dims          = 1:30,
    k             = 50,
    vote_min_frac = 0.55,             # majority needed to accept the label
    set_idents    = TRUE,             # write back to Idents()
    projected_prefix = "Projected: ", # final label prefix
    mouse_label_prefix_regex = "^Mouse:\\s*" # remove this from donor labels
) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Please install 'FNN'.")
  
  # pick first available reduction from choices
  reduction <- reduction[reduction %in% names(obj@reductions)][1]
  if (is.na(reduction)) stop("Chosen reduction not found in object.")
  
  emb_all <- Seurat::Embeddings(obj, reduction = reduction)
  
  # make sure dims exist
  dims <- dims[dims <= ncol(emb_all)]
  if (!length(dims)) stop("No valid dimensions found in the chosen reduction.")
  emb_all <- emb_all[, dims, drop = FALSE]
  
  # identify species
  stopifnot("species" %in% colnames(obj@meta.data))
  is_mouse <- obj$species == "Mouse"
  is_rat   <- obj$species == "Rat"
  if (!any(is_mouse) || !any(is_rat)) stop("No mouse or rat cells found via orig.ident.")
  
  mouse_cells <- which(is_mouse)
  rat_cells   <- which(is_rat)
  
  # donor labels (fine-grained) from current Idents
  mouse_labels <- as.character(Seurat::Idents(obj))[mouse_cells]
  
  # KNN from each rat cell to mouse cells
  nn <- FNN::get.knnx(
    data  = emb_all[mouse_cells, , drop = FALSE],
    query = emb_all[rat_cells,   , drop = FALSE],
    k     = k
  )$nn.index
  
  # majority vote per rat cell -> "Projected: XX"
  rat_proj <- vapply(seq_len(nrow(nn)), function(i) {
    labs <- mouse_labels[nn[i, ]]
    tab  <- sort(table(labs), decreasing = TRUE)
    top  <- names(tab)[1]
    frac <- tab[1] / sum(tab)
    if (frac >= vote_min_frac) {
      # strip any "Mouse: " (or your regex) then add "Projected: "
      top_clean <- sub(mouse_label_prefix_regex, "", top)
      paste0(projected_prefix, top_clean)
    } else NA_character_
  }, character(1))
  
  # store to meta
  out <- rep(NA_character_, ncol(obj)); names(out) <- colnames(obj)
  out[rat_cells] <- rat_proj
  obj$projected_label <- out
  
  # optionally set Idents: keep mouse as-is; replace rat where we have a projection
  if (set_idents) {
    cur <- as.character(Seurat::Idents(obj))
    cur[rat_cells] <- ifelse(is.na(rat_proj), cur[rat_cells], rat_proj)
    Seurat::Idents(obj) <- factor(cur)
  }
  
  obj
}

RatMouseTEC_EC.obj <- project_labels_from_mouse(
  RatMouseTEC_EC.obj,
  reduction     = "umap.cca",
  dims          = 1:50,
  k             = 1,
  vote_min_frac = 0.55,
  set_idents    = TRUE
)

DimPlot(subset(RatMouseTEC_EC.obj, subset =species =="Rat"), reduction = "umap.cca", label = T, label.box = T, repel = T)

table(Idents(subset(RatMouseTEC_EC.obj, subset =species =="Rat")))

v <- setNames(rep(NA_character_, ncol(RatMouseStroma.obj)),
              colnames(RatMouseStroma.obj))
v[names(RatMouseTEC_EC.obj$projected_label)] <- as.character(RatMouseTEC_EC.obj$projected_label)
RatMouseStroma.obj <- AddMetaData(RatMouseStroma.obj, metadata = v, col.name = "projected_label")

RatMouseStroma.obj$mixed_label <- RatMouseStroma.obj$projected_label

keep_ids <- c("Rat: capFb","Rat: TMC1","Rat: TMC2","Rat: TMC3","Rat: TMC4", "Rat: vSMC/PC")

id_all   <- Idents(RatMouseStroma.obj) 
id_chr   <- as.character(id_all)
cells_ok <- names(id_all)[id_all %in% keep_ids]

RatMouseStroma.obj$mixed_label[cells_ok] <- as.character(id_all[cells_ok])

RatMouseStroma.obj$mixed_label <- ifelse(is.na(RatMouseStroma.obj$mixed_label), 
                                         as.character(RatMouseStroma.obj$merged_cluster), RatMouseStroma.obj$mixed_label)

RatMouseStroma.obj$mixed_label <- factor(RatMouseStroma.obj$mixed_label,
                                         levels = c("Rat: capFb", "Rat: TMC1", "Rat: TMC2", "Rat: TMC3", "Rat: TMC4", "Rat: vSMC/PC",
                                                    "Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC","Projected: EC:fetal",
                                                    "Projected: MEC", "Projected: TEPC",
                                                    "Projected: TEC:early Pr", "Projected: TEC:cTEC", "Projected: TEC:mTEC1",
                                                    "Projected: TEC:mTEC-prol", "Projected: TEC:mTEC2",
                                                    "Projected: TEC:aaTEC1", "Projected: TEC:aaTEC2", 
                                                    "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic(basal)",
                                                    "Projected: TEC:mimetic(goblet)", "Projected: TEC:mimetic(muscle)",
                                                    "Projected: TEC:mimetic(neuroendo)", 
                                                    "Projected: TEC:mimetic(ionocyte)", "Projected: nmSC", 
                                                    "Mouse: FB:capsFB", "Mouse: FB:intFB", "Mouse: FB:medFB",  "Mouse: vSMC/PC", "Mouse: vSMC/PC:fetal", "Mouse: Fat", "Mouse: FB:fetal", 
                                                    "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC", "Mouse: MEC", "Mouse: EC:fetal",
                                                    "Mouse: TEPC", "Mouse: TEC:early Pr", "Mouse: TEC:cTEC", "Mouse: TEC:mTEC1", "Mouse: TEC:mTEC-prol", "Mouse: TEC:mTEC2",
                                                    "Mouse: TEC:aaTEC1", "Mouse: TEC:aaTEC2", 
                                                    "Mouse: TEC:mimetic(tuft)", "Mouse: TEC:mimetic(parathyroid)", "Mouse: TEC:mimetic(basal)", 
                                                    "Mouse: TEC:mimetic(microfold)", "Mouse: TEC:mimetic(goblet)", "Mouse: TEC:mimetic(muscle)", "Mouse: TEC:mimetic(neuroendo)",
                                                    "Mouse: TEC:mimetic(ciliated)", "Mouse: TEC:mimetic(ionocyte)", "Mouse: nmSC"))

                                                     
DimPlot(subset(RatMouseStroma.obj, subset=species %in% "Rat"), group.by = "mixed_label", 
        reduction = "umap.cca", label = T, label.box = T, repel = T, label.size = 2)+
  NoLegend() + ggtitle("Projected labels to Rat")
ggsave(filename = "pic/FigS2B_2.pdf", plot = get_last_plot(), width = 6, height = 5)

# Add projected clustering to the original rat thymocyte object
proj <- RatMouseTEC_EC.obj$projected_label
names(proj) <- colnames(RatMouseTEC_EC.obj)
proj <- proj[grep("^RatStroma_", names(proj))]
names(proj) <- sub("^RatStroma_", "", names(proj))
aligned <- proj[colnames(RatThymicStroma.obj)]
names(aligned) <- colnames(RatThymicStroma.obj)
RatThymicStroma.obj <- Seurat::AddMetaData(RatThymicStroma.obj, metadata = aligned, col.name = "projected_label")

DimPlot(RatThymicStroma.obj, group.by = "projected_label")

saveRDS(RatThymicStroma.obj, file = "20240331 rat thymic stroma scPCR recounted/RatThymicStroma.rds")

## remap endothelial cell subsets
DimPlot(subset(RatMouseStroma.obj, 
               subset = mixed_label %in% 
                 c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                   "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC")),
        reduction = "umap.cca", group.by = "species")+  xlim(-6.5, -2.5) + ylim(-14, -6)+
  ggtitle(NULL)+NoLegend()
ggsave(filename = "pic/Fig4A_1.pdf", plot = get_last_plot(), width = 5, height = 3)


DimPlot(subset(RatMouseStroma.obj, 
               subset = mixed_label %in% 
                 c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                   "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC")),
        reduction = "umap.cca", split.by = "species",
        group.by = "mixed_label", label = T, label.box = T)+
  xlim(-6.5, -2.5) + ylim(-14, -6)+ theme(strip.text = element_text(size = 20, face = "bold"))  + ggtitle(NULL)
ggsave(filename = "pic/Fig4A_2.pdf", plot = get_last_plot(), width = 12, height = 4)

FeaturePlot(subset(RatMouseStroma.obj, 
                   subset= mixed_label %in% 
                     c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                       "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC") ),
            reduction = "umap.cca", split.by  = "species", 
  features = c("Pecam1",
               "Selp","Vcam1","Icam1", # immigration TPEC
               "Bst1" # emigration TPE
               )) & coord_cartesian(xlim = c(-6.5, -2.5), ylim = c(-14, -6)) & theme(axis.title = element_blank())


DotPlot(subset(RatMouseStroma.obj, 
                   subset= mixed_label %in% 
                     c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                       "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC")),
        group.by = "mixed_label",
            features = rev(c("Pecam1",
                         "Selp","Vcam1","Icam1", # immigration TPEC
                         "Bst1" # emigration TPE
            )))+  
  coord_flip() +theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      legend.title = element_text(size = 10),
                      legend.text  = element_text(size = 9),
                      legend.key.height = unit(3, "mm"),
                      legend.key.width  = unit(3, "mm"),
                      legend.box.margin  = margin(t = 50, r = 0, b = 0, l = 0),
                      plot.margin        = margin(t = 2, r = 5, b = 5, l = 5)) +
  scale_y_discrete(limits = c("Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC",
                              "Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC"))
ggsave(filename = "pic/Fig4C.pdf", plot = get_last_plot(), width = 3.8, height = 2.9)

FeaturePlot(subset(RatMouseStroma.obj, 
                   subset= mixed_label %in% 
                     c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                       "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC") ),
            reduction = "umap.cca", split.by  = "species", 
            features = c("Sphk1", "Sphk2", "Spns2", "Sgpl1", "Plpp3" # S1p-associated
            )) & coord_cartesian(xlim = c(-6.5, -2.5), ylim = c(-14, -6)) & theme(axis.title = element_blank())

DotPlot(subset(RatMouseStroma.obj, 
               subset= mixed_label %in% 
                 c("Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC",
                   "Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC")),
        group.by = "mixed_label",
        features = rev(c("Sphk1", "Sphk2", "Spns2", "Sgpl1", "Plpp3" # S1p-associated
        )))+  
  coord_flip() +theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      legend.title = element_text(size = 10),
                      legend.text  = element_text(size = 9),
                      legend.key.height = unit(3, "mm"),
                      legend.key.width  = unit(3, "mm"),
                      legend.box.margin  = margin(t = 60, r = 0, b = 0, l = 0),
                      plot.margin        = margin(t = 2, r = 5, b = 5, l = 5)) +
  scale_y_discrete(limits = c("Mouse: EC:capEC", "Mouse: EC:vEC", "Mouse: EC:aEC",
                              "Projected: EC:capEC", "Projected: EC:vEC", "Projected: EC:aEC"))+
  scale_size_area(
    max_size = 5,                         # visual max dot (adjust)
    limits   = c(0, 100),                 # clamp range
    breaks   = c(1, 10, 100), labels = c("1", "10", "100"), trans="log1p")
ggsave(filename = "pic/Fig4D.pdf", plot = get_last_plot(), width = 3.8, height = 2.9)


#Ltbr controlled gene expression in Rat and Mouse stroma
PMID_32839611LtbrControlled <- intersect(readLines("PMID_32839611 thymic fibroblast/LtbrDelta_mFb vs flox_mFb.txt"),
                                         readLines("20240331 rat thymic stroma scPCR recounted/PMID_25224068 GNF Mouse GeneAtlas TRE list.txt")[-1])

PMID_32839611LtbrControlled_rat <- 
  as.character(convert_orthologs(as.data.frame(PMID_32839611LtbrControlled), input_species = "mouse", output_species = "rat", non121_strategy = 5, gene_input = "PMID_32839611LtbrControlled", gene_output = "dict"))

RatMouseStroma.obj <- AddModuleScore(object = RatMouseStroma.obj,
                                 features = list(PMID_32839611LtbrControlled_rat),
                                 name = "Ltbr controlled genes")
FeaturePlot(RatMouseStroma.obj,"Ltbr controlled genes1", reduction = "umap.cca", split.by = "species")
DotPlot(RatMouseStroma.obj,"Ltbr controlled genes1")

#remap mesenchymal cell subsets
unique(Idents(RatMouseStroma.obj))

RatMouseMesenchyme.obj <- subset(RatMouseStroma.obj, idents = c("Mouse: FB:capsFB", "Mouse: FB:intFB", "Mouse: FB:medFB", "Mouse: vSMC/PC", 
                                                                "Mouse: TEC:aaTEC2",
                                                        "Rat: TMC1", "Rat: TMC2", "Rat: TMC3", "Rat: capFb",
                                                        "Rat: TMC4", "Rat: vSMC/PC"))
DimPlot(RatMouseMesenchyme.obj, reduction = "umap.cca", split.by = "species")

DimPlot(RatMouseMesenchyme.obj, reduction = "umap.cca", group.by = "species")+ xlim(-13, 2) + ylim(-4, 6)+NoLegend()+
  ggtitle(NULL)
ggsave(filename = "pic/Fig3A_1.pdf", plot = get_last_plot(), width = 5, height = 4)

DimPlot(RatMouseMesenchyme.obj, reduction = "umap.cca", split.by = "species",
        label = T, label.box = T, repel = T)+ xlim(-13, 0) + ylim(-4, 6)+
  theme(strip.text = element_text(size = 20, face = "bold"))
ggsave(filename = "pic/Fig3A_2.pdf", plot = get_last_plot(), width = 12, height = 5)

# Fig. 3, 5.6 x 8 in
FeaturePlot(RatMouseMesenchyme.obj, features = c("Ltbr controlled genes1", "Cd34", "Pdgfra", "Pdgfrb"),
            reduction = "umap.cca", split.by = "species") & coord_cartesian(xlim = c(-13, 0), ylim = c(-4, 6)) & theme(axis.title = element_blank())
ggsave(filename = "pic/Fig3C.pdf", plot = get_last_plot(), width = 5.6, height = 8, device = cairo_pdf)

orth_MouseRat <- orthogene::report_orthologs(target_species = "mouse", reference_species = "rat", method_all_genes = "gprofiler", non121_strategy = "4") 

FeaturePlot(RatMouseMesenchyme.obj, features = c("Tslp", "Lgals5", "RT1-S3", "RT1-N3", "RT1-CE1"),  reduction = "umap.cca", split.by = "species") & # mouse mFb dominant
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
FeaturePlot(RatMouseMesenchyme.obj, features = c("Thy1", "Thbs1", "Thbs2", "Tgfb2", "Penk"),  reduction = "umap.cca", split.by = "species") & # Rat mFb dominant 1
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
FeaturePlot(RatMouseMesenchyme.obj, features = c("Nectin2", "Mdk", "Lpar1", "Lamb1", "Lama2"),  reduction = "umap.cca", split.by = "species") & # Rat mFb dominant 2
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
FeaturePlot(RatMouseMesenchyme.obj, features = c("Jag1", "Il6", "Il7", "Icoslg"),  reduction = "umap.cca", split.by = "species") & # Rat mFb dominant 3
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
FeaturePlot(RatMouseMesenchyme.obj, features = c("Fn1", "Cxcl12", "Ccl21", "Ccl19", "App", "Itgav", "Itgb1"),  reduction = "umap.cca", split.by = "species") & # Rat mFb dominant 4
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())

# Fig. 3, 4 x 7 in
DotPlot(subset(RatMouseStroma.obj, subset = mixed_label %in% c("Mouse: FB:medFB", "Mouse: TEC:mTEC2",
                                           "Rat: TMC3", "Rat: TMC4", "Projected: TEC:mTEC2")), 
        group.by = "mixed_label",
        features =  rev(c("Lgals5", "RT1-S3", "RT1-N3", "RT1-CE1", "Tslp",
                          "Ccl21", "Icoslg",
                          "Itgav", "Itgb1",
                          "App", "Cxcl12", "Fn1", "Lpar1", "Nectin2",
                          "Lama2", "Mdk", "Penk", "Tgfb2",
                          "Thbs1", "Il6", "Il7", "Jag1", "Lamb1", "Thbs2", "Thy1")))+
  RotatedAxis()+theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    coord_flip() +  scale_y_discrete(limits = c("Mouse: FB:medFB", "Mouse: TEC:mTEC2", "Rat: TMC3", "Rat: TMC4", "Projected: TEC:mTEC2")) +
  theme(
    legend.title = element_text(size = 12, angle = 90),
    legend.text  = element_text(size = 11, angle = 90)  # fallback for older ggplot2
  )
ggsave(filename = "pic/Fig3E.pdf", plot = get_last_plot(), width = 2.8, height = 7.8)

DotPlot(subset(Rat_MouseThymocyte.obj, idents=c("Projected: Mature CD4-1", "Mouse: Mature CD4-1")), 
        features =  c("Crlf2", "Ighm", "Ptprc", "Cd8b", "Adgre5", "Cd47", "Tgfbr1", "Tgfbr2", "Oprm1", "Cd226",
                      "Itga4", "Itga6", "Itgb1", "Notch1", "Il6r", "Il6st", "Il7r", "Il2rg", "Icos", "Cd28", 
                       "Itgb7", "Ccr7", "Sorl1"))+
  RotatedAxis()+theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Fig. 3, 5.6 x 5.7. in
FeaturePlot(RatMouseMesenchyme.obj, features = c("Lgals5", "RT1-CE1", "Tslp"),  reduction = "umap.cca", split.by = "species") & # mouse mFb dominant
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
ggsave(filename = "pic/Fig3F.pdf", plot = get_last_plot(), width = 5.6, height = 5.7)

# Fig. 3, 5.6 x 12 in
FeaturePlot(RatMouseMesenchyme.obj, features = c("Ccl21", "Icoslg", "Thbs1", "Il7"),  reduction = "umap.cca", split.by = "species") & # Rat mFb dominant
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
ggsave(filename = "pic/Fig3G.pdf", plot = get_last_plot(), width = 5.6, height = 7.8)

DotPlot(subset(RatMouseStroma.obj, subset=mixed_label %in% 
                 c("Mouse: FB:capsFB", "Mouse: FB:intFB", "Mouse: FB:medFB", "Mouse: vSMC/PC",
                   "Rat: capFb", "Rat: TMC1", "Rat: TMC2", "Rat: TMC3", "Rat: TMC4", "Rat: vSMC/PC")),
        group.by = "mixed_label",
        features =  rev(c("Sphk1","Sphk2", "Spns2", "Sgpl1", "Plpp3"
                      )), scale.max = NA)+
  RotatedAxis()+theme(axis.title = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
                      legend.title = element_text(size = 10),
                      legend.text  = element_text(size = 9),
                      legend.key.height = unit(3, "mm"),
                      legend.key.width  = unit(3, "mm"))+
  coord_flip() +  scale_y_discrete(limits = c("Mouse: FB:capsFB", "Mouse: FB:intFB", "Mouse: FB:medFB", "Mouse: vSMC/PC",
                                              "Rat: capFb", "Rat: TMC1", "Rat: TMC2", "Rat: TMC3", "Rat: TMC4", "Rat: vSMC/PC"))+
  scale_size_area(
    max_size = 5,                         # visual max dot (adjust)
    limits   = c(0.1, 100),                 # clamp range
    breaks   = c(0.1, 50, 100), labels = c("0.1", "50", "100"))
ggsave(filename = "pic/Fig3H.pdf", plot = get_last_plot(), width = 4.2, height = 3.2)

FeaturePlot(RatMouseMesenchyme.obj,
        split.by = "species",
        reduction = "umap.cca",
        features =  c("Sphk1","Sphk2", "Spns2", "Sgpl1", "Plpp3"))&
  coord_cartesian(xlim = c(-13, 2), ylim = c(-6, 6)) & theme(axis.title = element_blank())
ggsave(filename = "pic/Fig3I.pdf", plot = get_last_plot(), width = 5.6, height = 9.9)


# check rat or mouse only-genes caused by artifacts
mat <- GetAssayData(JoinLayers(RatMouseMesenchyme.obj), slot = "counts")
is_rat   <- RatMouseMesenchyme.obj$species == "Rat"
is_mouse <- RatMouseMesenchyme.obj$species == "Mouse"

det_rat   <- Matrix::rowSums(mat[, is_rat]   > 0) / sum(is_rat)
det_mouse <- Matrix::rowSums(mat[, is_mouse] > 0) / sum(is_mouse)

head(sort(det_rat[det_mouse == 0], decreasing = TRUE), 20)   # likely “rat-only” genes

features_common <- names(which(det_rat >0 & det_mouse >0)) # keep genes detected in ≥5% of cells in BOTH species (tune the 0.05)

# rTMC3 vs mmFb GSEA
TMC3_mFb_markers <- FindMarkers(JoinLayers(RatMouseMesenchyme.obj), ident.1 = c("Rat: TMC3"), ident.2 = c("Mouse: FB:medFB"),
                                only.pos=F, min.pct=0.1, logfc.threshold = 1,  features = features_common)

TMC3_mFb_markers <- TMC3_mFb_markers[order(TMC3_mFb_markers$avg_log2FC, decreasing=T),]

library(clusterProfiler)
library(enrichplot)
library(org.Rn.eg.db)

TMC3_mFb_list <- TMC3_mFb_markers$avg_log2FC
names(TMC3_mFb_list) <- rownames(TMC3_mFb_markers)

TMC3_mFb_list <- sort(TMC3_mFb_list, decreasing = T)

TMC3_mFb_GESA <- gseGO(geneList =TMC3_mFb_list, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL", pAdjustMethod="none")

TMC3_GESA <- TMC3_mFb_GESA; mFb_GESA <- TMC3_mFb_GESA
TMC3_GESA@result <- TMC3_GESA@result[TMC3_GESA@result$enrichmentScore > 0,]
mFb_GESA@result <- mFb_GESA@result[mFb_GESA@result$enrichmentScore < 0,]

treeplot(pairwise_termsim(TMC3_GESA)) | treeplot(pairwise_termsim(mFb_GESA))
set.seed(123)
emapplot(pairwise_termsim(TMC3_GESA), showCategory = 1200, group=T, nCluster = 5)

# rTMC4 vs mmFb GSEA
TMC4_mFb_markers <- FindMarkers(JoinLayers(RatMouseMesenchyme.obj), ident.1 = c("Rat: TMC4"), ident.2 = c("Mouse: FB:medFB"),
                                only.pos=F, min.pct=0.1, logfc.threshold = 1, features = features_common)

#TMC4_mFb_markers <- TMC4_mFb_markers[!(TMC4_mFb_markers$pct.1==0 | TMC4_mFb_markers$pct.2==0),]
TMC4_mFb_markers <- TMC4_mFb_markers[order(TMC4_mFb_markers$avg_log2FC, decreasing=T),]

TMC4_mFb_list <- TMC4_mFb_markers$avg_log2FC
names(TMC4_mFb_list) <- rownames(TMC4_mFb_markers)

TMC4_mFb_list <- sort(TMC4_mFb_list, decreasing = T)

TMC4_mFb_GESA <- gseGO(geneList =TMC4_mFb_list, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL"
                       #, pAdjustMethod="none"
                       )

dotplot(TMC4_mFb_GESA)

gseaplot2(TMC4_mFb_GESA, geneSetID = c("GO:0006939", "GO:0051707"),
          base_size = 15, subplots = 1:2, pvalue_table = T)

TMC4_GESA <- TMC4_mFb_GESA; mFb_GESA <- TMC4_mFb_GESA
TMC4_GESA@result <- TMC4_GESA@result[TMC4_GESA$enrichmentScore > 0,]
mFb_GESA@result <- mFb_GESA@result[mFb_GESA$enrichmentScore < 0,]

treeplot(pairwise_termsim(TMC4_GESA)) | treeplot(pairwise_termsim(mFb_GESA))

# rTMC3_4 vs mmFb GSEA
TMC3_4mFb_markers <- FindMarkers(JoinLayers(RatMouseMesenchyme.obj), ident.1 = c("Rat: TMC3", "Rat: TMC4"), ident.2 = c("Mouse: FB:medFB"),
                                only.pos=F, min.pct=0.1, logfc.threshold = 1, features = features_common)

TMC3_4mFb_markers <- TMC3_4mFb_markers[order(TMC3_4mFb_markers$avg_log2FC, decreasing=T),]

TMC3_4mFb_list <- TMC3_4mFb_markers$avg_log2FC
names(TMC3_4mFb_list) <- rownames(TMC3_4mFb_markers)

TMC3_4mFb_list  <- sort(TMC3_4mFb_list , decreasing = T)

TMC3_4mFb_GESA <- gseGO(geneList =TMC3_4mFb_list, OrgDb = "org.Rn.eg.db", ont = "ALL", keyType = "SYMBOL"
                        #, pAdjustMethod="none"
                        )

dotplot(TMC3_4mFb_GESA)

gseaplot2(TMC3_4mFb_GESA, geneSetID = c("GO:0006939", "GO:0051707"),
          base_size = 15, subplots = 1:2, pvalue_table = T)

TMC3_4GESA <- TMC3_4mFb_GESA; mFb_GESA <- TMC3_4mFb_GESA
TMC3_4GESA@result <- TMC3_4GESA@result[TMC3_4GESA$enrichmentScore > 0,]
mFb_GESA@result <- mFb_GESA@result[mFb_GESA$enrichmentScore < 0,]

emapplot(pairwise_termsim(TMC3_4GESA), showCategory=1200, group=T, nCluster = 5)

# GSEA about HALLMARK genes
library(msigdbr)

m_df <- msigdbr(species = "Mus musculus", category = "H")

m_df["r_gene_symbol"] <- convert_orthologs(m_df, input_species = "mouse", output_species = "rat", drop_nonorths = F,
                                           non121_strategy = 5, gene_input = "gene_symbol", gene_output = "column")[, "ortholog_gene"]

m_df <- rename(m_df, m_gene_symbol = gene_symbol)
m_df <- rename(m_df, gene_symbol = r_gene_symbol)

TMC3_4mFb_hallmark <- GSEA(geneList  = TMC3_4mFb_list,
                      TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
                      pvalueCutoff = 0.5)

dotplot(TMC3_4mFb_hallmark)

gseaplot2(TMC3_4mFb_hallmark, geneSetID = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
          base_size = 15, subplots = 1:2, pvalue_table = T)



# Integrate rat thymocyte and stroma objects----
RatThymicStroma.obj <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicStroma.rds")
RatThymocyte.obj <- readRDS("20250904 Rat thymocyte scRNA-seq/RatThymocyte.obj.rds")

RatThymocyte.obj <- RenameCells(RatThymocyte.obj, new.names = sub("-1$", "-2", colnames(RatThymocyte.obj)))

RatThymicRef <- merge(subset(RatThymicStroma.obj, ident= setdiff(levels(Idents(RatThymicStroma.obj)), c("Thymocyte", "B cell", "Myeloid"))),
                      subset(RatThymocyte.obj, ident= setdiff(levels(Idents(RatThymocyte.obj)), c("B+myeloid"))))

RatThymicRef <- RatThymicRef%>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

RatThymicRef <- IntegrateLayers(RatThymicRef, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
RatThymicRef <- FindNeighbors(RatThymicRef, reduction = "integrated.cca", dims = 1:30)
RatThymicRef  <- RunUMAP(RatThymicRef , reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(RatThymicRef, reduction = "umap.cca", label = T, label.box = T, repel = T)

RatThymicRef <- IntegrateLayers(RatThymicRef, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", assay = "RNA")
RatThymicRef <- FindNeighbors(RatThymicRef, reduction = "integrated.rpca", dims = 1:30,assay = "RNA")
RatThymicRef  <- RunUMAP(RatThymicRef , reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", assay = "RNA")
DimPlot(RatThymicRef, reduction = "umap.rpca", label = T, label.box = T, repel = T)

RatThymicRef <- IntegrateLayers(RatThymicRef, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
RatThymicRef <- FindNeighbors(RatThymicRef, reduction = "harmony", dims = 1:30)
RatThymicRef  <- RunUMAP(RatThymicRef , reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(RatThymicRef, reduction = "umap.harmony", label = T, label.box = T, repel = T)

RatThymicRef <- IntegrateLayers(RatThymicRef, method = JointPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.PCA")
RatThymicRef <- FindNeighbors(RatThymicRef, reduction = "integrated.PCA", dims = 1:30)
RatThymicRef  <- RunUMAP(RatThymicRef , reduction = "integrated.PCA", dims = 1:30, reduction.name = "umap.PCA")
DimPlot(RatThymicRef, reduction = "umap.PCA", label = T, label.box = T, repel = T)

RatThymicRef$consolidated_cluster <- RatThymicRef$projected_label

keep <- Idents(RatThymicRef) %in% c("capFb", "TMC1", "TMC2", "TMC3", "TMC4", "vSMC/PC")

RatThymicRef$consolidated_cluster[keep] <- as.character(Idents(RatThymicRef)[keep])

RatThymicRef$consolidated_cluster <- factor(RatThymicRef$consolidated_cluster, 
                                            levels = c("Projected: DN", "Projected: DP", "Projected: DP-P", 
                                                       "Projected: Pre DP-sig", "Projected: DP-sig", 
                                                       "Projected: Immature CD4-1", "Projected: Immature CD4-2", "Projected: intCD4", 
                                                       "Projected: Mature CD4-1", "Projected: Mature CD4-2", 
                                                       "Projected: Immature CD8-1", "Projected: Immature CD8-2", 
                                                       "Projected: Mature CD8-1", "Projected: Mature CD8-2", "Projected: Mature cycling T",
                                                       "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg",
                                                       "Projected: AgoSel2", "Projected: NKT", "Projected: gdT", 
                                                       "capFb", "TMC1", "TMC2", "TMC3", "TMC4", "vSMC/PC", 
                                                       "Projected: MEC", "Projected: nmSC", 
                                                       "Projected: TEPC", "Projected: TEC:early Pr", 
                                                       "Projected: TEC:cTEC", "Projected: TEC:mTEC1", 
                                                       "Projected: TEC:mTEC-prol", "Projected: TEC:mTEC2",
                                                       "Projected: TEC:aaTEC1", "Projected: TEC:aaTEC2",
                                                       "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic(goblet)", 
                                                       "Projected: TEC:mimetic(neuroendo)", "Projected: TEC:mimetic(muscle)", 
                                                       "Projected: TEC:mimetic(ionocyte)", "Projected: TEC:mimetic(basal)",
                                                       "Projected: EC:capEC", "Projected: EC:aEC", "Projected: EC:vEC", 
                                                       "Projected: EC:fetal", NA))

RatThymicRef <- subset(RatThymicRef, subset = !is.na(consolidated_cluster))

DimPlot(RatThymicRef, reduction = "umap.rpca", label = T, label.box = T, repel = T, group.by = "consolidated_cluster")+NoLegend()
DimPlot(subset(RatThymicRef, subset=consolidated_cluster %in% c( "Projected: AgoSel1", "Projected: Pre-Treg", "Projected: Treg",
                                                                 "Projected: AgoSel2", "Projected: NKT", "Projected: gdT", 
                                                                 "Projected: Immature CD4-1", "Projected: Immature CD4-2", "Projected: intCD4")),
        reduction = "umap.rpca", label = T, label.box = T, repel = T, group.by = "consolidated_cluster")+NoLegend()

saveRDS(RatThymicRef, "20240331 rat thymic stroma scPCR recounted/RatThymicRef.rds")
RatThymicRef <- readRDS("20240331 rat thymic stroma scPCR recounted/RatThymicRef.rds")

table(RatThymicRef$consolidated_cluster)

RatThymicRef_l <- subset(RatThymicRef, subset = !(consolidated_cluster %in%
                           c("Projected: intCD4", "Projected: Mature cycling T", "Projected: NKT", "Projected: gdT",
                             "Projected: MEC", "Projected: nmSC", "Projected: TEPC", "Projected: TEC:aaTEC1", "Projected: TEC:aaTEC2",
                             "Projected: TEC:mimetic(tuft)", "Projected: TEC:mimetic(goblet)", 
                             "Projected: TEC:mimetic(neuroendo)", "Projected: TEC:mimetic(muscle)", 
                             "Projected: TEC:mimetic(ionocyte)", "Projected: TEC:mimetic(basal)",
                             "Projected: EC:fetal")))

saveRDS(RatThymicRef_l, "20240331 rat thymic stroma scPCR recounted/RatThymicRef_l.rds")

DimPlot(RatThymicRef_l, reduction = "umap.rpca", label = T, label.box = T, repel = T, group.by = "consolidated_cluster")
ggsave(filename = "pic/FigS3_1.pdf", plot = get_last_plot(), width = 14, height = 7.5)

# Integrate Mouse thymocyte and stroma objects----
all_m_thymocyte <- readRDS("PMID_37580604 mouse thymocyte/all_m_thymocyte.rds")
MouseThymicStroma.obj <- readRDS(file = "PMID_39112630 mouse human integrated thymic stroma/MouseThymicStroma.rds")

WT_m_thymocyte <- subset(all_m_thymocyte,  
                                   subset= genotype %in% "WT") 

WT_m_thymocyte@meta.data$Thymo_Stroma <- "Thymocyte"

Ad_MouseThymicStroma.obj <- subset(MouseThymicStroma.obj,  
                               subset=
                                 #age_range %in% "2-14wk" & # remove fetal & old subsets 
                                genotype %in% c("wild type")) 

Ad_MouseThymicStroma.obj@meta.data$Thymo_Stroma <- "Stroma"

MouseThymicAll.obj <- merge(WT_m_thymocyte, Ad_MouseThymicStroma.obj)

MouseThymicAll.obj <- MouseThymicAll.obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

MouseThymicAll.obj <- IntegrateLayers(MouseThymicAll.obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
MouseThymicAll.obj <- FindNeighbors(MouseThymicAll.obj, reduction = "integrated.cca", dims = 1:50)
MouseThymicAll.obj  <- RunUMAP(MouseThymicAll.obj, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")

DimPlot(MouseThymicAll.obj, reduction = "umap.cca", label = T, label.box = T)+NoLegend()

saveRDS(MouseThymicAll.obj, "PMID_39112630 mouse human integrated thymic stroma/MouseThymicAll.rds")

# cell-cell interaction analysis (CellChat)----
library(CellChat)
# Rat
cellchat <- readRDS("20240331 rat thymic stroma scPCR recounted/cellchat.rds")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

pdf("incoming_heatmap.pdf", width = 7, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 7)
dev.off()

pdf("outgoing_heatmap.pdf", width = 7, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 2)
dev.off()

lr_out <- subsetCommunication(cellchat, sources.use = c("TMC4"), targets.use = "Projected: Mature CD4-1")

# Fig. 3, 3.2 x 9.2 in
netVisual_bubble(
  cellchat,
  sources.use   = c("TMC3","TMC4", "Projected: TEC:mTEC1", "Projected: TEC:mTEC2"),
  targets.use   = c("Projected: Mature CD4-1", "Projected: Mature CD4-2"),
  remove.isolate = TRUE,
  angle.x        = 90,
  font.size      = 12)
ggsave(filename = "pic/Fig3D_2.pdf", plot = get_last_plot(), width = 3.2, height = 9.2)

netVisual_bubble(
  cellchat,
  sources.use   = c("TMC3","TMC4", "Projected: TEC:mTEC1", "Projected: TEC:mTEC2"),
  targets.use   = c("Projected: AgoSel1"),
  remove.isolate = TRUE,
  angle.x        = 90,
  font.size      = 12)
ggsave(filename = "pic/Fig3D_2.pdf", plot = get_last_plot(), width = 3.2, height = 9.2)

netVisual_bubble(
  cellchat,
  sources.use   = c("Projected: Mature CD4-1", "Projected: Mature CD4-2"),
  targets.use   = c("Projected: TEC:mTEC2","TMC3","TMC4"),
  remove.isolate = TRUE,
  angle.x        = 45,
  font.size      = 11)

# Mouse
cellchat_m <- readRDS("PMID_39112630 mouse human integrated thymic stroma/cellchat_m.rds")

cellchat_m@idents <- factor(cellchat_m@idents,
                            levels = c("DN", "DP-P", "DP", "Pre DP-sig", 
                                       "DP-sig", "Immature CD4-1", "Immature CD4-2", "intCD4", "Mature CD4-1", 
                                       "Mature CD4-2", "Immature CD8-1", "Immature CD8-2", "Mature CD8-1", 
                                       "Mature CD8-2", "Mature cycling T", "AgoSel1", "Pre Treg", "Treg", 
                                       "AgoSel2", "NKT", "gdT", "Erythrocyte", "Myeloid", "B cell", "MEC",
                                       "Fat", "FB:capsFB", "FB:intFB", "FB:medFB", "FB:fetal", "vSMC/PC", "vSMC/PC:fetal", 
                                       "TEPC",  "TEC:early Pr", "TEC:cTEC", "TEC:mTEC1", "TEC:mTEC-prol", "TEC:mTEC2", 
                                       "TEC:aaTEC1", "TEC:aaTEC2", 
                                       "TEC:mimetic(tuft)", "TEC:mimetic(basal)", "TEC:mimetic(goblet)", "TEC:mimetic(microfold)", 
                                       "TEC:mimetic(ciliated)", "TEC:mimetic(parathyroid)", "TEC:mimetic(muscle)", "TEC:mimetic(neuroendo)", "nmSC", 
                                       "EC:capEC", "EC:aEC","EC:vEC", "EC:fetal"))


cellchat_m <- computeCommunProbPathway(cellchat_m)

cellchat_m <- aggregateNet(cellchat_m)

groupSize <- as.numeric(table(cellchat_m@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_m@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_m@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_m@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_m, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_m, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat_m, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_m, signaling = pathways.show, color.heatmap = "Reds")

cellchat_m <- netAnalysis_computeCentrality(cellchat_m, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat_m, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat_m, pattern = "incoming")

pdf("incoming_heatmap.pdf", width = 7, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_m, pattern = "incoming", font.size = 7)
dev.off()

pdf("outgoing_heatmap.pdf", width = 7, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_m, pattern = "outgoing", font.size = 2)
dev.off()

lr_out <- subsetCommunication(cellchat_m, sources.use = c("FB: medFB"), targets.use = "Mature CD4-1")

# Fig. 3, 3.2 x 6.5 in
netVisual_bubble(
  cellchat_m,
  sources.use   = c("FB:medFB"),
  targets.use   = "Mature CD4-1",
  remove.isolate = TRUE,
  font.size      = 15)
ggsave(filename = "pic/Fig3D_1.pdf", plot = get_last_plot(), width = 3.2, height = 6.5)

netVisual_bubble(
  cellchat_m,
  sources.use   = c("FB:medFB", "TEC:aaTEC2", "TEC:mTEC2"),
  targets.use   = c("Mature CD4-1", "Mature CD4-2", "Mature CD8-1", "Mature CD8-2"),
  remove.isolate = TRUE,
  font.size      = 15)


#make combined table
FbmTECvsCD4_m <- CellChat::subsetCommunication(cellchat_m, 
                              sources.use = c("FB:medFB","TEC:aaTEC2", "TEC:mTEC2"), 
                              targets.use = c("Mature CD4-1", "Mature CD4-2"))

FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Ccl21a  - Ccr7"] <- "Ccl21  - Ccr7"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Lgals9  - Cd45"] <- "Lgals5  - Cd45"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Lgals9  - Ighm"] <- "Lgals5  - Ighm"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Icosl  - Icos"] <- "Icoslg  - Icos"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Icosl  - Icos"] <- "Icoslg  - Icos"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Icosl  - Cd28"] <- "Icoslg  - Cd28"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="Jam1  - (Itgal+Itgb2)"] <- "F11r  - (Itgal+Itgb2)"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-m2  - Cd8b1"] <- "RT1-M2  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-m3  - Cd8b1"] <- "RT1-M3-1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-q10  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-q7  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-t22  - Cd8b1"] <- "RT1-N3  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-t23  - Cd8b1"] <- "RT1-S3  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-q6  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-d1  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-k1  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-q4  - Cd8b1"] <- "RT1-CE1  - Cd8b"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-aa  - Cd4"] <- "RT1-Ba  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-ab1  - Cd4"] <- "RT1-Bb  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-eb1  - Cd4"] <- "RT1-Db1  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-dma  - Cd4"] <- "RT1-DMa  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-dmb1  - Cd4"] <- "RT1-DMb  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-dmb2  - Cd4"] <- "RT1-DMb  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-oa  - Cd4"] <- "RT1-DOa  - Cd4"
FbmTECvsCD4_m$interaction_name_2[FbmTECvsCD4_m$interaction_name_2 =="H2-ob  - Cd4"] <- "RT1-DOb  - Cd4"

FbmTECvsCD4_r <- CellChat::subsetCommunication(cellchat,
                                   sources.use = c("TMC3","TMC4", "Projected: TEC:mTEC2"),
                                   targets.use = c("Projected: Mature CD4-1", "Projected: Mature CD4-2"))

FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="Icosl  - Cd28"] <- "Icoslg  - Cd28"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="Icosl  - Icos"] <- "Icoslg  - Icos"
FbmTECvsCD4_r <- FbmTECvsCD4_r[!FbmTECvsCD4_r$interaction_name_2 =="Ccl21  - Ccr7",]
FbmTECvsCD4_r <- FbmTECvsCD4_r[!FbmTECvsCD4_r$interaction_name_2 =="Ccl21b  - Ccr7",]
FbmTECvsCD4_r <- FbmTECvsCD4_r[!FbmTECvsCD4_r$interaction_name_2 =="Ccl21d  - Ccr7",]
FbmTECvsCD4_r <- FbmTECvsCD4_r[!FbmTECvsCD4_r$interaction_name_2 =="Ccl21f  - Ccr7",]
FbmTECvsCD4_r <- FbmTECvsCD4_r[!FbmTECvsCD4_r$interaction_name_2 =="Ccl21e  - Ccr7",]
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="Ccl21a  - Ccr7"] <- "Ccl21  - Ccr7"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="Jam1  - (Itgal+Itgb2)"] <- "F11r  - (Itgal+Itgb2)"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-m3  - Cd8a"] <- "RT1-M3-1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q1  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q10  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q7  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-t23  - Cd8a"] <- "RT1-S3  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q6  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-d1  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-k1  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q9  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-T26  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-q4  - Cd8a"] <- "RT1-CE1  - Cd8a"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-m3  - Cd8b1"] <- "RT1-M3-1  - Cd8b"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-ea-ps  - Cd4"] <- "RT1-Da  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-aa  - Cd4"] <- "RT1-Ba  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-ab1  - Cd4"] <- "RT1-Bb  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-eb1  - Cd4"] <- "RT1-Db1  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-dma  - Cd4"] <- "RT1-DMa  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-dmb1  - Cd4"] <- "RT1-DMb  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-dmb2  - Cd4"] <- "RT1-DMb  - Cd4"
FbmTECvsCD4_r$interaction_name_2[FbmTECvsCD4_r$interaction_name_2 =="H2-oa  - Cd4"] <- "RT1-DOa  - Cd4"

library(dplyr)

norm_comm <- function(df) {
  # pick/construct pair name
  pair_col <-
    if ("interaction_name_2" %in% names(df)) "interaction_name_2" else
      if ("interaction_name"   %in% names(df)) "interaction_name"   else
        if (all(c("ligand","receptor") %in% names(df))) {
          df$interaction_name_2 <- paste(df$ligand, df$receptor, sep = " - ")
          "interaction_name_2"
        } else stop("No interaction name column (or ligand/receptor) in this table.")
  
  # probability column (name differs across versions)
  prob_col <-
    if ("prob" %in% names(df)) "prob" else
      if ("prob.comb" %in% names(df)) "prob.comb" else
        stop("No probability column found (looked for 'prob' or 'prob.comb').")
  
  # p-value column (many variants)
  pcol <- intersect(c("pval","p.value","p","pvalue"), names(df))
  if (length(pcol) == 0) df$pval <- NA_real_ else df$pval <- df[[pcol[1]]]
  
  # standardize + add the x-axis label used in bubble plots
  df %>%
    transmute(
      source, target,
      pair   = .data[[pair_col]],
      prob   = .data[[prob_col]],
      pval,
      x_label = paste(source, "->", target)
    )
}

d_m <- norm_comm(FbmTECvsCD4_m)
d_r <- norm_comm(FbmTECvsCD4_r)

combined_long <- bind_rows(d_m, d_r)

library(ggplot2)

# order of columns on x-axis
cols <- c("FB:medFB -> Mature CD4-1", "FB:medFB -> Mature CD4-2",
          "TEC:mTEC2 -> Mature CD4-1", "TEC:mTEC2 -> Mature CD4-2",
          
          "TMC3 -> Projected: Mature CD4-1","TMC3 -> Projected: Mature CD4-2",
          "TMC4 -> Projected: Mature CD4-1", "TMC4 -> Projected: Mature CD4-2",
          "Projected: TEC:mTEC2 -> Projected: Mature CD4-1", "Projected: TEC:mTEC2 -> Projected: Mature CD4-2",
          "TEC:aaTEC2 -> Mature CD4-1", "TEC:aaTEC2 -> Mature CD4-2")

d_plot <- combined_long %>%
  mutate(x_label = factor(x_label, levels = cols),
         p_show = ifelse(is.finite(-log10(pval)), -log10(pval), NA_real_),
         shape_flag = ifelse(!is.na(pval) & pval < 0.01, "p<0.01", "ns"))%>%
  mutate(pval = suppressWarnings(as.numeric(pval)),
         # categorize p-values; treat 0 as p<0.01
         p_cat = case_when(
           !is.na(pval) & pval <= 0.01 ~ "p<0.01",
           !is.na(pval) & pval <= 0.05 ~ "0.01<p<=0.05",
           TRUE                                ~ "ns"
         ), .keep = "unused")

# order of columns on Y-axis
pair_order <- c("Ccl19  - Ccr7", "Icam1  - (Itgal+Itgb2)", "Icam1  - Itgal", 
                "Mdk  - Ncl", "Ptn  - Ncl", "RT1-Bb  - Cd4", "Icam1  - Spn", "Lgals5  - Cd45", 
                "Lgals5  - Ighm", "RT1-CE1  - Cd8b", "RT1-N3  - Cd8b", 
                "RT1-S3  - Cd8b", "Tslp  - (Il7r+Crlf2)", "Ccl21  - Ccr7", 
                "Icoslg  - Cd28", "Icoslg  - Icos", "(Itgav+Itgb1) - Adgre5", 
                "App - Sorl1",  
                "Fn1  - (Itga4+Itgb7)", "Lpar1 - Adgre5", "Nectin2  - Cd226", 
                "RT1-M3-1  - Cd8a", "RT1-CE1  - Cd8a", "Lrrc4c - Ptprf",
               "RT1-S3  - Cd8a", "Mdk  - (Itga4+Itgb1)", "Fn1  - (Itga4+Itgb1)", "Mdk  - (Itga6+Itgb1)", "Penk  - Oprm1", "Tgfb2  - (Tgfbr1+Tgfbr2)", 
               "Lamb2  - (Itga6+Itgb1)", "Lama2  - (Itga6+Itgb1)", 
               
               "Igf1  - Igf1r", "Thbs1  - Cd47", "Cxcl12  - Cxcr4", 
                "Il6  - (Il6r+Il6st)", "Il7  - (Il7r+Il2rg)", "Jag1  - Notch1", 
                "Thbs2  - Cd47", 
                "Thy1 - Adgre5", "Lamb1  - (Itga6+Itgb1)", 
                "Grn  - Sort1", 
                "Cd70  - Cd27", "Cd80  - Cd28", "F11r  - (Itgal+Itgb2)", 
                "Nectin1  - Cd96", "RT1-Ba  - Cd4",  
                "RT1-Db1  - Cd4", "RT1-DMa  - Cd4", "RT1-DMb  - Cd4", 
                "RT1-DOa  - Cd4", "RT1-M3-1  - Cd8b",  
                "Gdf15  - Tgfbr2", "RT1-DOb  - Cd4", "RT1-M2  - Cd8b", "Ccl25  - Ccr9",
                "Cd80  - Cd274", "Cd86  - Cd28", "Lama3  - (Itga6+Itgb1)", 
                "Lamc2  - (Itga6+Itgb1)", "Pvr  - Cd226", "RT1-Da  - Cd4", "Tgfb1 - (Tgfbr1+Tgfbr2)")

d_plot <- d_plot %>% mutate(pair = as.character(pair)) %>%
  mutate(pair = factor(pair, levels = rev(c(pair_order))))

ggplot(d_plot, aes(x = x_label, y = pair)) +
  geom_point(aes(size = p_cat, fill = prob, shape = p_cat),
             stroke = 0.3, colour = "black") +
  scale_fill_gradientn(name = "Commun. Prob.",
                       colours = c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"),
                       na.value = "white") +
  scale_size_manual(name = "p-value", values = size_map, drop = FALSE) +
  scale_shape_manual(name = "p-value", values = shape_map, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank())
ggsave(filename = "pic/Fig3D.pdf", plot = get_last_plot(), width = 4, height = 12)

