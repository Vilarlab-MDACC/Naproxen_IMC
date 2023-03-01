library(tidyverse)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(Seurat)
library(ggpubr)
library(M3C)
library(ComplexHeatmap)
library(circlize)

read_csv("meta_sample.csv")%>%
  right_join(read_csv("merge.csv"), by =c("...1" = "Metadata_roi...6")) -> cell_data


stain_with_mask = cell_data[25:61]%>%rename_with(~str_remove_all(.x,"Intensity_MeanIntensity_"))

colon <- CreateSeuratObject(counts = stain_with_mask%>%t%>%`colnames<-`(1:ncol(.)), project = "IMC")
all.genes <- rownames(colon)
colon <- NormalizeData(colon,normalization.method = "CLR", margin = 1)
colon <- ScaleData(colon, features = all.genes)


out<-Rphenograph(stain_with_mask)
clusters <- factor(membership(out[[2]]))
colon$pg_clusters <- clusters



GetAssayData(colon)%>%
  as_tibble(rownames = NA)%>%
  rownames_to_column("feature")%>%
  pivot_longer(-1,
               names_to="cell",
               values_to = "density")%>%
  left_join(colon$pg_clusters%>%
              as_tibble_col("cluster")%>%
              rownames_to_column("cell"),
            by="cell")%>%
  group_by(feature,cluster)%>%
  summarise(
    mean = mean(density),
    median = median(density)
  )%>%
  ungroup%>%
  group_by(feature)%>%
  mutate(
    mean_scaled = scale(mean)[,1] ,
    median_scaled = scale(median)[,1]
  ) -> cluster_result



col_order = c("CD3e","CD45","CD19","CD4","CD8a","GranzymeB","OX40","CD137","ICOS","Ki67","CD45RO","CCR6","CD11b","CD11c","CD14","CD66b","CD68","TIM3","VISTA","B7-H3","PD-L1","PD-1", "GAPDH","CK","NaK-ATPase","HLA-DR","B2M","CD39","CD73","MPO","CD163","FoxP3")
row_order = c(149, 124, 130, 24)
col_annotation = factor(c(rep("Lymphocytes Makers",5),rep("T Cell Activation Makers",5),rep("T Cell Memory Makers",2),rep("Monocyte Markers",5),rep("Immune Check Points",5),rep("Epithelial Makers",3),rep("Other",7)),
                        c("Lymphocytes Makers","T Cell Activation Makers","T cell Memory Makers","Monocyte Markers","Immune Check Points","Epithelial Markes","Other"))


cluster_result%>%
  filter(cluster%in%c(149, 124, 130, 24))%>%
  group_by(feature)%>%
  pivot_wider(id_cols = cluster,
              names_from = "feature",
              values_from = "mean_scaled")%>%
  column_to_rownames("cluster")%>%
  select(col_order)   %>%mutate_all(~replace(.,.>1,1))-> df.to_plt


df.to_plt %>%
  Heatmap(
    column_order = col_order,
    row_order = c("149", "124", "130", "24"),
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, f) {
      grid.circle(
        x = x,
        y = y,
        r = min(unit.c(width, height)) * 2,
        gp = gpar(fill = col_fun(.[i, j]), col = NA)
      )
    },
    col = function(x)
      colorRamp2(c(-0.5, 0.5), c("yellow", "blue"))(x),
    bottom_annotation = HeatmapAnnotation(`Marker Category` = col_annotation),
    heatmap_legend_param = list(title = "Expression"),
    show_row_names =F ,
    border = "black",
    left_annotation =  rowAnnotation(population = anno_text(
      c(
        "CD11c+ IDO-1+ PD-L1+ Macrophage",
        "HLA-DR+ CD8+ T cells",
        "CD11c- Macrophages",
        "Proliferating Effector CD8+ T cells"
        
      ),
      just = "right"    ,location = 1)),
    right_annotation = rowAnnotation(
      Compartment = anno_text_box(c("b","b","b","a"),
                                  list(a="Epitheium",
                                       b="Out of Epitheium"),
                                  by = "anno_block",
                                  side = "right"
                                  
      )
      
    ))
