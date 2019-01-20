# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)
}

# Usage example
packages<-c("optparse", "circlize")
check.packages(packages)

#install.packages('optparse')
#install.packages('circlize')
#suppressWarnings(suppressMessages(library(optparse)))
#suppressWarnings(suppressMessages(library(circlize)))

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

png(filename=opt$plot, units="in", width=10, height=10, pointsize=12, res=300)
grid.col = c(A = 'brown1', B = 'lawngreen', C = 'mediumorchid', D = 'mediumslateblue', E = 'royalblue', F = 'sandybrown')
par(mar = rep(0,4), cex = 1.2)
chordDiagram(t(mat), grid.col = grid.col)
#circos.track(circos.text, facing = "clockwise")
invisible(dev.off())

