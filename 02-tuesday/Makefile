all: slides birds ponds dune pdf

slides: slides.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("slides.Rmd")'

ponds: ponds.Rmd
	Rscript -e "rmarkdown::render(\"ponds.Rmd\")"

birds: birds.Rmd
	Rscript -e "rmarkdown::render(\"birds.Rmd\")"

dune: dune-meadows.Rmd
	Rscript -e "rmarkdown::render(\"dune-meadows.Rmd\")"

pdf: slides.Rmd slides.html
	Rscript -e "renderthis::to_pdf(from = \"slides.Rmd\")"
