PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check clean

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	Rscript -e 'pkgdown::build_site()'
	Rscript -e 'pkgdown::build_articles()'
	Rscript -e 'rmarkdown::render("README.Rmd")'

readme2:
	Rscript -e 'rmarkdown::render("README.Rmd", "html_document")'

build: rd
	R CMD build $(PKGSRC)

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: 
	Rscript -e 'devtools::check()'

check2: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	$(RM) -r $(PKGNAME).Rcheck/

