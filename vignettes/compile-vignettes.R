library(knitr)

rmarkdown::render("vignette-basics.Rmd", output_file="vignette-basics.html")
rmarkdown::render("vignette-basics.Rmd", output_file="vignette-basics.pdf")

rmarkdown::render("vignette-ploidy.Rmd", output_file="vignette-ploidy.html")
rmarkdown::render("vignette-ploidy.Rmd", output_file="vignette-ploidy.pdf")

rmarkdown::render("vignette-rCNAs.Rmd", output_file="vignette-rCNAs.html")
rmarkdown::render("vignette-rCNAs.Rmd", output_file="vignette-rCNAs.pdf")
