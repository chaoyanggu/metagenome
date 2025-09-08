library("ggplot2")

# 获取当前目录下所有以.length结尾的文件
length_files <- list.files(pattern = "\\.length.GC$")

# 遍历每个文件
for (file in length_files) {
  # 读取文件内容到数据框
  data <- read.table(file, header = FALSE)
  colnames(data) <- c("contig","length","GC")
  
  # 检查第二列是否存在
  if (ncol(data) < 2) {
    stop(paste("File", file, "does not contain enough columns."))
  }
  
  # 为当前文件创建标题
  plot_title <- paste("Density Plot for", file)
  
  data1 <- data[data$contig!="Total",]
  # 绘制density图
  p <- ggplot(data1, aes(x = length)) + geom_density(fill = "blue", alpha = 0.5) + labs(title = plot_title, x = "Length", y = "Density") +theme_classic()
  
  # 打印图形到控制台
  print(p)
  
  # 保存图形为PNG文件
  ggsave(paste0(file, ".density.pdf"), plot = p, width = 10, height = 6)
}

