all: slides purl ponds meadows permanova dbrda baxter pdf

slides: slides.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("slides.Rmd")'

purl: slides.Rmd
	Rscript -e "knitr::purl(\"slides.Rmd\")"

ponds: constrained-ordination.Rmd
	Rscript -e "rmarkdown::render(\"constrained-ordination.Rmd\")"

meadows: spring-meadows.Rmd
	Rscript -e "rmarkdown::render(\"spring-meadows.Rmd\")"

permanova: permanova.Rmd
	Rscript -e "rmarkdown::render(\"permanova.Rmd\")"

dbrda: dbrda.Rmd
	Rscript -e "rmarkdown::render(\"dbrda.Rmd\")"

baxter: baxter-otus.Rmd
	Rscript -e "rmarkdown::render(\"baxter-otus.Rmd\")"

pdf: slides.Rmd slides.html
	Rscript -e "renderthis::to_pdf(from = \"slides.Rmd\")"
