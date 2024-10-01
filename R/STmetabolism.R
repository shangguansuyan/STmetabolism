#' Title ST.metabolism Function
#' @import Seurat
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom scMetabolism sc.metabolism.Seurat
#' @importFrom dplyr %>%
#' @importFrom CellTrek traint celltrek
#' @importFrom Seurat RenameCells Cells
#'
#' @title Single-cell example data (sc)
#' @description A small example dataset of single-cell RNA-seq data.
#' @format A Seurat object with 25314 features and 8944 samples.
"sc"
#'
#' @title Spatial example data (st)
#' @description A small example dataset of spatial transcriptomics data.
#' @format A Seurat object with 36601 features and 2261 samples.
"st"
#'
#' @title ST.metabolism Function
#' @description This function integrates single-cell and spatial data for
#'     metabolic analysis.
#' @param st_data 一个V4空间转录组Seurat对象
#' @param sc_data 一个V4单细胞转录组Seurat对象
#' @param sc_assay 单细胞的assay对象名称
#' @param cell_names metadata中细胞名称的列名
#' @param int_assay 嵌入空间转录组数据的单细胞数据的assay对象
#' @param reduction 选择降维的方式
#' @param intp 如果为TRUE则执行插值
#' @param intp_pnt 插值点数
#' @param intp_lin 如果为TRUE则执行线性插值
#' @param nPCs PC数量
#' @param ntree 树数
#' @param dist_thresh 阈值距离
#' @param top_spot 一个单元格可以绘制的最大点数
#' @param spot_n 一个点可以包含的最大单元格数
#' @param repel_r 排斥半径
#' @param repel_iter 排斥迭代
#' @param keep_model 如果为TRUE则返回经过训练的随机森林模型
#' @param method 进行单细胞代谢评分的方法
#' @param imputation 是否在代谢评分之前插补其数据
#' @param ncores 计算单细胞代谢使用的内核数
#' @param metabolism.type 选择代谢分析的参考数据库
#'
#' @return 一个列表，包含STmetbolism对象，celltrek对象和训练后的嵌入数据集
#' @export
#' @name ST.metabolism
#' @examples
#' # library(Seurat,lib.loc = 'D:\Program Files\R\R-4.4.1\seurat_lib')
#'     #V4版本的Seurat
#' library(CellTrek)
#' library(scMetabolism)
#' library(dplyr)
#' # library(ggplot2)
#' # library(rsvd)
#' library(STmetabolism)
#' # 需要输入V4 Seurat创建的单细胞对象和空转对象
#' head(sc)
#' sc
#' head(st)
#' st
#' ST_met <- ST.metabolism(st_data = st, sc_data = sc)
#' # 查看训练后的嵌入数据
#' DimPlot(ST_met$traint, group.by = "type")
#' # celltrek对象可以直接按照celltrek的教程运行后续步骤
#' CellTrek::celltrek_vis(
#'   ST_met$celltrek@meta.data %>%
#'     dplyr::select(coord_x, coord_y, cell_type:id_new),
#'   ST_met$celltrek@images$anterior1@image,
#'   ST_met$celltrek@images$anterior1@scale.factors$lowres
#' )
#' # 以celltrek的点为基础绘制切片的代谢活性分布图
#' SpatialFeaturePlot(ST_met$STmetbolism_obj,
#'   features = c("Glycolysis / Gluconeogenesis"), pt.size.factor = 2.5
#' )
#'
ST.metabolism <- function(st_data,
                          sc_data,
                          sc_assay = "RNA",
                          cell_names = "cell_type",
                          int_assay = "traint",
                          reduction = "pca",
                          intp = TRUE,
                          intp_pnt = 5000,
                          intp_lin = FALSE,
                          nPCs = 30,
                          ntree = 1000,
                          dist_thresh = 0.55,
                          top_spot = 5,
                          spot_n = 5,
                          repel_r = 20,
                          repel_iter = 20,
                          keep_model = TRUE,
                          method = "VISION",
                          imputation = FALSE,
                          ncores = 1,
                          metabolism.type = "KEGG") {
  # 重命名细胞名称
  st_data <- RenameCells(st_data, new.names = make.names(Cells(st_data)))
  sc_data <- RenameCells(sc_data, new.names = make.names(Cells(sc_data)))

  # 运行 CellTrek traint 和 celltrek
  traint <- CellTrek::traint(
    st_data = st_data, sc_data = sc_data,
    sc_assay = sc_assay, cell_names = cell_names
  )
  celltrek <- CellTrek::celltrek(
    st_sc_int = traint, int_assay = int_assay,
    sc_data = sc_data, sc_assay = sc_assay,
    reduction = reduction, intp = intp,
    intp_pnt = intp_pnt, intp_lin = intp_lin,
    nPCs = nPCs, ntree = ntree,
    dist_thresh = dist_thresh,
    top_spot = top_spot,
    spot_n = spot_n, repel_r = repel_r,
    repel_iter = repel_iter,
    keep_model = keep_model
  )$celltrek

  # 对整合后的celltrek对象创建Seurat对象
  countexp.Seurat <- CreateSeuratObject(
    counts = celltrek@assays$RNA@data,
    meta.data = celltrek@meta.data
  )

  # 进行代谢分析
  countexp.Seurat <- sc.metabolism.Seurat(
    obj = countexp.Seurat,
    method = method,
    imputation = imputation,
    ncores = ncores,
    metabolism.type = metabolism.type
  )

  # 提取代谢矩阵
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  meta_metabolism <- as.data.frame(t(metabolism.matrix))

  # 创建ST对象并整合代谢数据
  ST <- sc_data
  ST@assays[["Spatial"]] <- celltrek@assays[["RNA"]]
  ST@meta.data <- meta_metabolism
  ST@assays[["RNA"]] <- celltrek@assays[["RNA"]]
  ST@active.ident <- celltrek@active.ident
  ST@graphs <- celltrek@graphs
  ST@images <- celltrek@images

  # 创建最后的返回对象
  STmet <- list(
    STmetbolism_obj = ST, celltrek_obj = celltrek,
    traint_obj = traint
  )

  # 清理内存
  rm(
    countexp.Seurat, meta_metabolism,
    metabolism.matrix, ST, traint, celltrek
  )
  gc()

  # 返回结果
  return(STmet)
}
