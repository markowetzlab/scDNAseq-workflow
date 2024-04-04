library(knitr)

rmarkdown::render("vignette-basics.Rmd", output_file="vignette-basics.html")
rmarkdown::render("vignette-basics.Rmd", output_file="vignette-basics.pdf")

rmarkdown::render("vignette-ploidy", output_file="vignette-ploidy")
rmarkdown::render("vignette-ploidy", output_file="vignette-ploidy")

rmarkdown::render("vignette-rCNAs.Rmd", output_file="vignette-rCNAs.html")
rmarkdown::render("vignette-rCNAs.Rmd", output_file="vignette-rCNAs.pdf")
