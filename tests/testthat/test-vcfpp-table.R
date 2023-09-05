library(vcfppR)

vcffile <- system.file("extdata", "test-GL.vcf.gz", package="vcfppR")
(gt <- tableGT(vcffile,"chr20"))
(gl <- tableGL(vcffile,"chr20"))

vcffile <- system.file("extdata", "test-PL.vcf.gz", package="vcfppR")
(pl <- tablePL(vcffile,"chr20"))

