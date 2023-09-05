#include "vcfpp.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' parse GT FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of genotypes for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableGT(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getRegionIndex(region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
    IntegerVector pos(nsnps);
    NumericVector qual(nsnps);
    vector<vector<int>> GT(nsnps);
    vector<int> gt;
    for (int i = 0; i < nsnps; i++) {
        vcf.getNextVariant(var);
        pos(i) = var.POS();
        qual(i) = var.QUAL();
        chr(i) = var.CHROM();
        id(i) = var.ID();
        ref(i) = var.REF();
        alt(i) = var.ALT();
        filter(i) = var.FILTER();
        info(i) = var.INFO();
        var.getGenotypes(gt);
        GT[i] = gt;
    }
    return List::create(Named("samples") = vcf.header.getSamples(),
                        Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                        Named("ref") = ref, Named("alt") = alt,
                        Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

//' parse GP FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of genotype posterior probabilites for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableGP(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<float>> GP(nsnps);
  vector<float> gp;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("GP", gp);
    GP[i] = gp;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("gp") = GP);
}

//' parse DS FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of genotype dosages for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableDS(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<float>> DS(nsnps);
  vector<float> ds;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("DS", ds);
    DS[i] = ds;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("ds") = DS);
}

//' parse GL FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of log-scaled genotype likelihoods for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableGL(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<float>> GL(nsnps);
  vector<float> gl;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("GL", gl);
    GL[i] = gl;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("gl") = GL);
}

//' parse PL FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of phred-scaled genotype likelihoods for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tablePL(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> PL(nsnps);
  vector<int> pl;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("PL", pl);
    PL[i] = pl;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("pl") = PL);
}

//' parse GQ FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of phred-scaled genotype quality for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableGQ(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> GQ(nsnps);
  vector<int> gq;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("GQ", gq);
    GQ[i] = gq;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("gq") = GQ);
}


//' parse HQ FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of phred-scaled haplotype quality for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableHQ(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> HQ(nsnps);
  vector<int> hq;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("HQ", hq);
    HQ[i] = hq;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("hq") = HQ);
}

//' parse DP FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of read depths for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableDP(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> DP(nsnps);
  vector<int> dp;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("DP", dp);
    DP[i] = dp;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("dp") = DP);
}

//' parse MQ FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of RMS mapping quality for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tableMQ(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> MQ(nsnps);
  vector<int> mq;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("MQ", mq);
    MQ[i] = mq;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("mq") = MQ);
}

//' parse PQ FORMAT of the VCF file into tables in R
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of phasing quality for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List tablePQ(std::string vcffile, std::string region,
             std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region,samples);
  vcfpp::BcfRecord var(vcf.header);
  int nsnps = vcf.getRegionIndex(region);
  CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps),
      info(nsnps);
  IntegerVector pos(nsnps);
  NumericVector qual(nsnps);
  vector<vector<int>> PQ(nsnps);
  vector<int> pq;
  for (int i = 0; i < nsnps; i++) {
    vcf.getNextVariant(var);
    pos(i) = var.POS();
    qual(i) = var.QUAL();
    chr(i) = var.CHROM();
    id(i) = var.ID();
    ref(i) = var.REF();
    alt(i) = var.ALT();
    filter(i) = var.FILTER();
    info(i) = var.INFO();
    var.getFORMAT("PQ", pq);
    PQ[i] = pq;
  }
  return List::create(Named("samples") = vcf.header.getSamples(),
                      Named("chr") = chr, Named("pos") = pos, Named("id") = id,
                      Named("ref") = ref, Named("alt") = alt,
                      Named("qual") = qual, Named("filter") = filter,
                      Named("info") = info, Named("pq") = PQ);
}
