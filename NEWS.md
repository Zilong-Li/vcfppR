# vcfppR 0.8.3
* bugfix: vcfcomp with truth being an vcftable object for stats = "r2"!

# vcfppR 0.8.2
* extend subset S3 function for vcftable
* vcfcomp: stats='gtgq'
* fix vcfplot
* more unit tests and examples

# vcfppR 0.8.1
* remove unused code in vcfpp.h
* vcfcomp: rename the col name "concordance" as "accuracy"

# vcfppR 0.8.0
* copy static libs into libs folder

# vcfppR 0.7.6
* fix CRAN notes

# vcfppR 0.7.5
* fix issue #14 to handle libdeflate dependency on linux

# vcfppR 0.7.4
* fix issue #13 to speedup setRegion

# vcfppR 0.7.3
* enable libdeflate on linux

# vcfppR 0.7.2
* make `htslib` available for including and linking by other packages

# vcfppR 0.7.0
* fix parsing BCF when subsetting sample ids and the genomic region at the same time
* fix `vcfcomp` with `stats="nrc"` when the input vectors have different length
* use the latest version of htslib v1.21
* reduce the size of installation

# vcfppR 0.6.0
* upgrade `vcfpp.h` to v0.6.0
* new API `vcfreader$setRegion`
* new API `vcfreader$getStatus`

# vcfppR 0.5.0
* upgrade `vcfpp.h` to v0.5.1
* add `bcfreader$updateSamples()` function
* `vcfcomp`: can take `vcftable` object as input
* `vcfcomp`: defaults changed `setid=FALSE`
* bug fix


# vcfppR 0.4.6
* add `vcfplot` function
* fix bug in `vcftable` when FORMAT has unaligned list 
* add vignette on `vcfcomp`

# vcfppR 0.4.5

* add `vcfcomp` function
* fix memory leaks and pass ASAN check

# vcfppR 0.4.0

* add `setid` option for vcftable
* patches for upcoming Rtools on windows

# vcfppR 0.3.8

* fix issues on M1 Mac

# vcfppR 0.3.7

* add `vcfreader@rmFormatTag()`
* add `vcfreader@samples()`
* `vcftable` supports vartype = 'sv'

# vcfppR 0.3.6

* add `vcfreader::line`
* remove `vcfreader::setVariant`
* bug fix `setFormatStr`
* more units tests

# vcfppR 0.3.5

* add vcfreader and vcfwriter

# vcfppR 0.3.4

* support filters by variant ID

# vcfppR 0.3.3

* use API of vcfpp.h v0.3.1
* set missing vaules in FORMAT to NA
* support filters by QUAL, FILTER and INFO

# vcfppR 0.3.0

* First release using vcfpp.h v0.3.0
* add vcftable function
* add vcfsummary function
* add vcfpopgen function
