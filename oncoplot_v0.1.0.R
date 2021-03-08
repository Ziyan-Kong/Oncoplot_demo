setwd("D:/JobManagement/工作任务/20210304001-基因组突变景观图/")


# 加载软件包
library(ComplexHeatmap)
library(openxlsx)
library(magrittr)

# Demo 数据准备
mut <- read.xlsx("mutation.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE) %>% t %>% as.data.frame
rownames(mut) <- lapply(rownames(mut), function(x){
                          gene <- strsplit(x, split = ".", fixed = TRUE)[[1]][1]
                        }) %>% unlist
mut[is.na(mut)] <- ""

head(mut[1:5, 1:3])


# 突变类型颜色标记
col <- c("nonsynonymous SNV"="#FF0000",
         "frameshift insertion"="#FF00FF", 
         "frameshift deletion"="#0000FF", 
         "stopgain"="#00FFFF", 
         "nonframeshift deletion"="#008000",
         "splicing"="#00FA9A", 
         "nonframeshift substitution"="#D2B48C"
)
## the characterize of mutation
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "nonsynonymous SNV" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["nonsynonymous SNV"], col = NA))
  },
  "frameshift insertion" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["frameshift insertion"], col = NA))
  },
  "frameshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["frameshift deletion"], col = NA))
  },
  "stopgain" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["stopgain"], col = NA))
  },
  "splicing" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["splicing"], col = NA))
  },
  "nonframeshift substitution" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["nonframeshift substitution"], col = NA))
  },
  "nonframeshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["nonframeshift deletion"], col = NA))
  }
)

# clinical data
phe <- read.xlsx("clinical.xlsx", sheet = 1, rowNames = FALSE, colNames = TRUE, check.names = FALSE, sep.names = " ")

# patients to clinical characteristics
phe <- subset(phe, phe$`sample ID` %in% colnames(mut))
phe$`mRNA subtype` <- ifelse(phe$`mRNA subtype`==1, "SI", 
                               ifelse(phe$`mRNA subtype`==2, "SII", "SIII"))
phe$`Proteomic subtype` <- ifelse(phe$`Proteomic subtype`==1, "SI", 
                               ifelse(phe$`Proteomic subtype`==2, "SII", "SIII"))
phe$`mRNA subtype` <- factor(x = phe$`mRNA subtype`, levels = sort(unique(phe$`mRNA subtype`)))
phe$`Proteomic subtype` <- factor(x = phe$`Proteomic subtype`, levels = sort(unique(phe$`Proteomic subtype`)))

phe[1:5, 1:5]

mut <- mut[, phe$`sample ID`]



# the type of mutation
hlp <- list(
  title = "Alternations", 
  at = c("nonsynonymous SNV","frameshift insertion", "frameshift deletion",
         "stopgain", "nonframeshift deletion","splicing", "nonframeshift substitution"), 
  labels = c("nonsynonymous SNV","frameshift insertion", "frameshift deletion",
             "stopgain", "nonframeshift deletion","splicing", "nonframeshift substitution"),
  direction = "horizontal",
  nrow = 3
)

# annotation for clinical characteristics
han <- HeatmapAnnotation(
  "cbar" = anno_oncoprint_barplot(),
  "Age" = phe$Age,
  "Gender" = phe$Gender,
  "mRNA Subtype" = phe$`mRNA subtype`,
  "Proteomic Subtype" = phe$`Proteomic subtype`,
  "Tumor Purity" = phe$`Tumor purity`,
  "TNM Stage" = phe$`TNM stage`,
  "OS" = phe$`Overall survial`,
  col = list(
    "mRNA Subtype" = c("SI"="orange", "SII"="green", "SIII"="skyblue"),
    "Proteomic Subtype" = c("SI"="orange", "SII"="green", "SIII"="skyblue"),
    "Gender" = c("female"="red", "male"="blue")
  ),
  annotation_name_side = "left"
)

anno <- oncoPrint(mut, 
          alter_fun = alter_fun, 
          col = col, 
          heatmap_legend_param = hlp,
          remove_empty_columns = TRUE,
          remove_empty_rows = TRUE,
          row_names_side = "left",
          pct_side = "right",
          top_annotation = han
  )
draw(anno, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")

