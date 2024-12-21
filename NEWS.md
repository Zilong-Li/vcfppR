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
