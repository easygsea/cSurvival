library(circlize)
library(RColorBrewer)

display.brewer.all()

categories <- c("DepMap","TARGET","TCGA")
categories_cols = brewer.pal(3, "Set2")

projectss <- unname(unlist(projects[categories]))

circos.par("track.height" = 0.3)
circos.initialize(categories, xlim = cbind(rep(0, 3), sapply(categories,function(x) length(projects[[x]]))))
circos.track(categories, ylim = c(0, 1), track.index = 1, track.height = 0.3, panel.fun = function(x, y) {
  circos.axis(h = 1, major.tick = F, minor.ticks = F,
              labels.cex = 0.1, col = categories_cols[CELL_META$sector.numeric.index], labels.col="#ffffff")
  circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, 
              facing = "bending.inside",col = categories_cols[CELL_META$sector.numeric.index])
}, bg.border = NA)
circos.clear()

# par(new = TRUE) # <- magic
# circos.par("canvas.xlim" = c(-2, 1), "canvas.ylim" = c(-2, 1))
# circos.initialize(projectss, xlim = c(0,1))
# circos.track(projectss, ylim = c(-1, 1), track.index = 1, track.height = 0.1)
# circos.clear()

