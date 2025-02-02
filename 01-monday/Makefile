all: slides tatry vltava cluster

slides: slides.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("slides.Rmd")'

tatry: tatry.Rmd
	Rscript -e "rmarkdown::render(\"tatry.Rmd\")"

vltava: vltava.Rmd
	Rscript -e "rmarkdown::render(\"vltava.Rmd\")"

cluster: cluster-analysis.Rmd
	Rscript -e "rmarkdown::render(\"cluster-analysis.Rmd\")"