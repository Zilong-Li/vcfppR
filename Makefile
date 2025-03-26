PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check clean

rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme: rd
	Rscript -e 'rmarkdown::render("README.Rmd")'
	Rscript -e 'pkgdown::build_site(examples=FALSE)'
	Rscript -e 'pkgdown::build_articles()'

readme2: rd
	Rscript -e 'rmarkdown::render("README.Rmd", "html_document")'

build: rd
	cd ..;\
	R CMD build $(PKGSRC)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check:
	cd ..;\
	# Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz", args="--as-cran")'
	Rscript -e 'devtools::check()'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz  --as-cran

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

