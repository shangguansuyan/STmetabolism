# STmetabolism
整合celltrek和scMetabolism结果以单细胞精度映射至切片上

library(devtools)
install_github("navinlabcode/CellTrek")
devtools::install_github("wu-yc/scMetabolism")
install_github("shangguansuyan/STmetabolism")

library(Seurat) #V4版本的Seurat和SeuratObject
library(CellTrek)
library(scMetabolism)
library(dplyr)
library(ggplot2)
library(rsvd)
library(STmetabolism)

# 使用示例数据
sc<-sc #V4版本的单细胞Seurat对象
st<-st #V4版本的空转Seurat对象
ST_met <- ST.metabolism(st_data = st, sc_data = sc)

# 查看训练后的嵌入数据
DimPlot(ST_met$traint, group.by = "type")

# celltrek对象可以直接按照celltrek的教程运行后续步骤
CellTrek::celltrek_vis(
  ST_met$celltrek@meta.data %>%
    dplyr::select(coord_x, coord_y, cell_type:id_new),
  ST_met$celltrek@images$anterior1@image,
  ST_met$celltrek@images$anterior1@scale.factors$lowres
)

# 以celltrek的点为基础绘制切片的代谢活性分布图
SpatialFeaturePlot(ST_met$STmetbolism_obj,
  features = c("Glycolysis / Gluconeogenesis"), pt.size.factor = 2.5
)

