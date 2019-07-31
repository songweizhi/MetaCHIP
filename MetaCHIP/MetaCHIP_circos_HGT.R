
# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)
}


# install packages if not 
packages<-c("optparse", "circlize")
invisible(suppressMessages(check.packages(packages)))


options(warn=-1)

option_list = list(
  
  make_option(c("-m", "--matrix"), 
              type="character", 
              help="input matrix", 
              metavar="character"),
  
  make_option(c("-p", "--plot"), 
              type="character", 
              help="output plot", 
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mat = read.table(opt$matrix, header = TRUE)

png(filename=opt$plot, units="in", width=25, height=25, pointsize=12, res=300)
grid.col = c(A = 'brown1', B = 'lawngreen', C = 'mediumorchid', D = 'mediumslateblue', E = 'royalblue', F = 'sandybrown')
par(mar = rep(0,4), cex = 1.2)
label_order = sort(union(rownames(mat), colnames(mat)))
chordDiagram(t(mat), order = label_order, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)

# rorate label
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

#circos.track(circos.text, facing = "clockwise")
invisible(dev.off())

rm(list=ls())
