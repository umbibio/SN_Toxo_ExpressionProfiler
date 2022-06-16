library(shiny)
library(DT)
library(Seurat)
library(ggplot2)
library(rsconnect)


#S.O.RH348 <- readRDS('S_O_RH348.rds')
S.O.Toxo <- readRDS('S.O.intra_lables.rds')
S.O.SN <- readRDS('S.O_sn_toxo_labs.rds')
S.O.PB <- readRDS('S.O.pb.labs.rds')
genes <- readRDS('toxo_sn_pb_prod_desc.rds')


S.O.Toxo@meta.data$phase <- factor(S.O.Toxo@meta.data$phase, 
                                   levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
Idents(S.O.Toxo) <- 'phase'

S.O.SN@meta.data$seurat_clusters <- factor(S.O.SN@meta.data$seurat_clusters,
                                           levels = sort(unique(S.O.SN@meta.data$seurat_clusters)))

Idents(S.O.SN) <- 'seurat_clusters'

S.O.PB@meta.data$cells <- factor(S.O.PB@meta.data$cells, 
                                 levels = c("RngE", "RngL", "TrpE", "TrpM", "TrpL", "SchE", "SchM", "SchL"))

Idents(S.O.PB) <- 'cells'

server <- function(input, output, session) {

  df <- as.data.frame(genes)
  options(DT.options = list(pageLength = 5))
  
  # row selection
  output$x7 = DT::renderDataTable(df, server = FALSE, selection = 'single')
  
  ## expression  
  output$x1 = renderPlot({
    s = input$x7_rows_selected
    selected.GT1 <- gsub('_', '-', df$Tg_gene_id[s])
    if(length(s)){
      if(!is.na(selected.GT1)){
        if(selected.GT1 %in% rownames(S.O.Toxo)){
          p1 = FeaturePlot(object = S.O.Toxo, features = selected.GT1, 
                           label = T, pt.size = 1, label.size = 8,
                           cols = c("lightgrey", "blue"), reduction = "umap") + 
            theme_minimal() + 
            labs(title="T. gondii",x ="UMAP1", y = "UMAP2") + 
            theme(
              plot.title = element_text(size=14, face="bold.italic"),
              axis.title.x = element_text(size=12, face="bold"),
              axis.title.y = element_text(size=12, face="bold")
            )
          
          plot(p1)
        }
      }
    }
  })
  
  output$x2 = renderPlot({
    s = input$x7_rows_selected
    selected.SN <- gsub('_', '-', df$sn_gene_id[s])
    if(length(s)){
      if(!is.na(selected.SN)){
        if(selected.SN %in% rownames(S.O.SN)){
          p1 = FeaturePlot(object = S.O.SN, features = selected.SN, 
                           label = T, pt.size = 1, label.size = 8,
                           cols = c("lightgrey", "blue"), reduction = "umap") + 
            theme_minimal() + 
            labs(title="S. Neurona",x ="UMAP1", y = "UMAP2") + 
            theme(
              plot.title = element_text(size=14, face="bold.italic"),
              axis.title.x = element_text(size=12, face="bold"),
              axis.title.y = element_text(size=12, face="bold")
            )
          plot(p1)
        }
      }
    }
  })
  
  output$x3 = renderPlot({
    s = input$x7_rows_selected
    selected.PB <- gsub('_', '-', df$pb_gene_id[s])
    if(length(s)){
      if(!is.na(selected.PB)){
        if(selected.PB %in% rownames(S.O.PB)){
          p1 = FeaturePlot(object = S.O.PB, features = selected.PB, 
                           label = T, pt.size = 1, label.size = 8,
                           cols = c("lightgrey", "blue"), reduction = "umap") + 
            theme_minimal() + 
            labs(title="P. Berghei ",x ="UMAP1", y = "UMAP2") + 
            theme(
              plot.title = element_text(size=14, face="bold.italic"),
              axis.title.x = element_text(size=12, face="bold"),
              axis.title.y = element_text(size=12, face="bold")
            )
          plot(p1)
        }
      }
    }
  })
  
  
  ## Vln 
  output$x4 = renderPlot({
    s = input$x7_rows_selected
    selected.GT1 <- gsub('_', '-', df$Tg_gene_id[s])
    if(length(s)){
      if(!is.na(selected.GT1)){
        if(selected.GT1 %in% rownames(S.O.Toxo)){
          p1 = VlnPlot(object = S.O.Toxo, features = selected.GT1) + 
            theme_minimal() + 
            labs(title="T. gondii",x ="", y = "log2(expr)") + 
            NoLegend() + 
            theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(
              plot.title = element_text(size=14, face="bold"),
              axis.title.x = element_text(size=14, face="bold", hjust = 1),
              axis.title.y = element_text(size=14, face="bold")
            )
          
          
          plot(p1)
        }
      }
    }
  })
  
  output$x5 = renderPlot({
    s = input$x7_rows_selected
    selected.SN <- gsub('_', '-', df$sn_gene_id[s])
    if(length(s)){
      if(!is.na(selected.SN)){
        if(selected.SN %in% rownames(S.O.SN)){
          p1 = VlnPlot(object = S.O.SN, features = selected.SN) + 
            theme_minimal() + 
            labs(title="S. Neurona",x ="", y = "log2(expr)") + 
            NoLegend() + 
            theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(
              plot.title = element_text(size=14, face="bold"),
              axis.title.x = element_text(size=14, face="bold", hjust = 1),
              axis.title.y = element_text(size=14, face="bold")
            )
          
          plot(p1)
        }
      }
    }
  })
  
  output$x6 = renderPlot({
    s = input$x7_rows_selected
    selected.PB <- gsub('_', '-', df$pb_gene_id[s])
    if(length(s)){
      if(!is.na(selected.PB)){
        if(selected.PB %in% rownames(S.O.PB)){
          p1 = VlnPlot(object = S.O.PB, features = selected.PB) + 
            theme_minimal() + 
            labs(title="P. Berghei ",x ="", y = "log2(expr)") + NoLegend() + 
            theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
            theme(
              plot.title = element_text(size=14, face="bold"),
              axis.title.x = element_text(size=14, face="bold", hjust = 1),
              axis.title.y = element_text(size=14, face="bold")
            )
          
          plot(p1)
        }
      }
    }
  })
  
  session$allowReconnect(TRUE)
}