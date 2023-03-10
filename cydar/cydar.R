library(tidyverse)
library(edgeR)
library(cydar)
library(ncdfFlow)
library(flowCore)
library(ggpubr)
library(ComplexHeatmap)

markers_loc = c(4:38)

list.files("Cells/","csv",full.names = T)%>%
  map(
    ~ read_csv(.x)%>%
      dplyr::filter(maskGreen>=128)%>%
      select(markers_loc%>%all_of)%>%
      list%>%
      `names<-`(.x%>%
                  str_remove(".csv")%>%
                  str_remove("Cells//"))%>%
      poolCells
  )%>%
  `names<-`(list.files("Cells/","csv",full.names = T)%>%
              str_remove(".csv")%>%
              str_remove("Cells//"))%>%
  as("flowSet")%>%
  ncdfFlowSet -> collected.exprs.ep



list.files("Cells/","csv",full.names = T)%>%
  map(
    ~ read_csv(.x)%>%
      dplyr::filter(maskGreen<128)%>%
      select(markers_loc%>%all_of)%>%
      list%>%
      `names<-`(.x%>%
                  str_remove(".csv")%>%
                  str_remove("Cells//"))%>%
      poolCells
  )%>%
  `names<-`(list.files("Cells/","csv",full.names = T)%>%
              str_remove(".csv")%>%
              str_remove("Cells//"))%>%
  as("flowSet")%>%
  ncdfFlowSet -> collected.exprs.nep


cydar.process <- function(collected.exprs) {
  pool.ff <- poolCells(collected.exprs)
  trans <- estimateLogicle(pool.ff, colnames(pool.ff))
  processed.exprs <- transform(collected.exprs, trans)
  cd <- prepareCellData(processed.exprs)
  cd <- countCells(cd)
  return(cd)
}

cd.ep <- collected.exprs.ep %>% cydar.process()
cd.nep <- collected.exprs.nep %>% cydar.process()




edgeR.process <- function(cd){
  y <- DGEList(assay(cd), lib.size = cd$totals)
  keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
  cd <- cd[keep, ]
  y <- y[keep, ]
  
  colnames(y) %>%
    as_tibble %>%
    mutate(
      case = str_extract(value, "case_\\d+"),
      naproxen = str_extract(value, "Naproxen_(High|Low|Placebo)"),
      visit = str_extract(value, "(Pre|Post)-treatment"),
      roi = str_extract(value, "roi_\\d*")
    ) -> meta
  
  design <- model.matrix( ~ 0+factor(meta$case))
  naproxen_post = meta$naproxen != "Naproxen_Placebo" &
    meta$visit == "Post-treatment"
  placebo_post = meta$naproxen == "Naproxen_Placebo" &
    meta$visit == "Post-treatment"
  
  
  design <- cbind(design, naproxen_post, placebo_post)
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, contrast = c(rep(0, ncol(design) - 2), -1, 1))
  qvals <- spatialFDR(intensities(cd), res$table$PValue)
  
  return(list(cd = cd ,res =res,qvals = qvals))
}
 
res.ep = cd.ep%>%edgeR.process() 
res.ep %>% saveRDS("res.ep.rds")
res.nep = cd.nep%>%edgeR.process()
res.nep %>% saveRDS("res.nep.rds")

 





wilcox.process <- function(cd) {
  colnames(cd) %>%
    as_tibble %>%
    mutate(
      case = str_extract(value, "case_\\d+"),
      naproxen = str_extract(value, "Naproxen_(High|Low|Placebo)"),
      visit = str_extract(value, "(Pre|Post)-treatment"),
      roi = str_extract(value, "roi_\\d*")
    ) -> meta
  
  
    assay(cd) %>%
      as_tibble(rownames = NA) %>%
      rownames_to_column("sphere") %>%
      pivot_longer(starts_with("case")) %>%
      left_join(meta, by = c("name" = "value")) %>%
      mutate(group = factor(
        case_when(str_detect(naproxen, "Placebo") ~ "Placebo",
                  T ~ "Naproxen"),
        levels = c("Placebo", "Naproxen")
      )) %>%
      group_by(case, group, visit, roi) %>%
      mutate(proportion = value / sum(value)) %>%
      group_by(case, sphere, group, visit) %>%
      summarise(proportion = mean(proportion)) %>%
      pivot_wider(
        id_cols = c(case, sphere, group),
        values_from = proportion,
        names_from = visit
      ) %>%
      mutate(diff = `Post-treatment` - `Pre-treatment`,
             rel_diff = diff / (`Pre-treatment`+1e-6)) %>%
      group_by(sphere) %>%
      group_map(~ as_tibble_row(
        list(
          sphere = .y$sphere,
          p.rel_diff = wilcox.test(rel_diff ~ group, data = .x)$p.value,
          p.diff = wilcox.test(diff ~ group, data = .x)$p.value
        )
      )) %>%
      bind_rows() %>%
      mutate(
        fdr.diff = p.adjust(p.diff, method = "fdr"),
        fdr.rel_diff = p.adjust(p.rel_diff, method = "fdr"),
      ) %>% return
  }
  
  
wilcox.ep = cd.ep %>% wilcox.process 
wilcox.ep%>%write_rds("wilcox.ep.rds")

wilcox.nep = cd.nep %>% wilcox.process 
wilcox.ep%>%write_rds("wilcox.nep.rds")
