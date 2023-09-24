#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List tableGT(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

// [[Rcpp::export]]
List tableGP(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gp") = GP);
}

// [[Rcpp::export]]
List tableDS(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("ds") = DS);
}

// [[Rcpp::export]]
List tableGL(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gl") = GL);
}

// [[Rcpp::export]]
List tablePL(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("pl") = PL);
}

// [[Rcpp::export]]
List tableGQ(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gq") = GQ);
}

// [[Rcpp::export]]
List tableHQ(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("hq") = HQ);
}

// [[Rcpp::export]]
List tableDP(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("dp") = DP);
}

// [[Rcpp::export]]
List tableMQ(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("mq") = MQ);
}

// [[Rcpp::export]]
List tablePQ(std::string vcffile, std::string region, std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsnps = vcf.getVariantsCount(var, region);
    CharacterVector chr(nsnps), ref(nsnps), alt(nsnps), id(nsnps), filter(nsnps), info(nsnps);
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
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("pq") = PQ);
}
