# PSO_GMM_scRNA_clustering
## High dimention(rely on your dataset) scRNA dataset cell clustering tool (naive edition)

### 加载包
library(PSOClustering)

### 步骤1：加载数据
data <- load_data("path/to/fibro.rds", "path/to/deg.xlsx")

### 步骤2：初始化PSO参数
pso_params <- initialize_pso(data$deg, size = 30, generation = 3, toprange = 15)

### 步骤3：运行PSO算法
best_deg_solution <- run_pso(pso_params$deg_list, size = 30, generation = 3)

### 步骤4：运行GMM聚类
gmm_result <- run_gmm(data$fibro, n_clusters = 3)

### 查看结果
print(gmm_result$clusters)
