
rm(list = ls())
# globals
path.site <- "C:/Users/Devin/Dropbox/Projects/dincerti.github.io"

# convert r markdown to markdown
rmd2md <- function(path_site = path.site,
                   file, dir_rmd = "_rmd", dir_figs = "figs/",
                   out_ext='.md', in_ext='.rmd',
                   recursive = FALSE) {
  require(knitr)
  content <- readLines(file.path(path_site, dir_rmd, file))
  outfile <- file.path(path_site, paste0(substr(file, 1, (nchar(file)-(nchar(in_ext)))), out_ext))
  render_jekyll(highlight = "pygments")
  opts_knit$set(out.format='markdown') 
  opts_knit$set(base.url = "/", base.dir = path_site, root.dir = path_site)
  opts_chunk$set(fig.path = dir_figs, fig.width = 8, fig.height = 5, fig.align='center')
  knit(text = content, output = outfile)
}

# list of files
files <- list.files(path = file.path(path.site, "_rmd"), pattern = ".rmd", 
                    ignore.case = TRUE, recursive = FALSE)

# convert
rmd2md(file = "twopart.Rmd")
rmd2md(file = "dynamic_twopart.Rmd")
rmd2md(file = "diabetes_highcost.Rmd")
rmd2md(file = "markov_cohort.Rmd")
rmd2md(file = "bayesian_markov_cohort.Rmd")
rmd2md(file = "bayesian_meta_analysis.Rmd")
