# Rscript saw_sunburst.R [values]

args<-commandArgs(T)

library(plotly)
library(dplyr)

f <- list(
  family = "Arial")

fig <- plot_ly(
  labels = c("Total Reads", "Invalid CID Reads", "Valid CID Reads", "Non-Relevant Short Reads", "Clean Reads", "Unmapped Reads", "Mapped Reads", "Multi-Mapped Reads", "Uniquely-Mapped Reads", "Unannotated Reads", "Annotated Reads", "Unique Reads", "Duplicated Reads"),
  parents = c("", "Total Reads", "Total Reads", "Valid CID Reads", "Valid CID Reads", "Clean Reads", "Clean Reads", "Mapped Reads", "Mapped Reads", "Uniquely-Mapped Reads", "Uniquely-Mapped Reads", "Annotated Reads", "Annotated Reads"),
  values = c(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13]),
  type = 'sunburst',
  marker = list(colors = list("red", "yellow", "orange", "palegreen", "green", "rosybrown", "salmon", "skyblue", "royalblue", "thistle", "violet", "cyan", "gray")),
  branchvalues = 'total'
) %>% layout(font=f)

#export(fig, file = "stats_sunburst.png")

htmlwidgets::saveWidget(as_widget(fig), paste0("Sequencing","_sunburst.html"), title=paste0("Plot", " Key Sequencing Metrices"))
