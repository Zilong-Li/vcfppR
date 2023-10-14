#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List tableGT(std::string vcffile, std::string region, std::string samples = "-", double qualval = 0,
             bool pass = false, bool INFO = true, bool snps = false, bool indels = false,
             bool multiallelics = false, bool multisnps = false) {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<vector<int>> GT;
    vector<int> gt, pos;
    vector<float> qual;
    vector<std::string> chr, ref, alt, id, filter, info;
    while (vcf.getNextVariant(var)) {
        if (pass && (var.FILTER() != "PASS")) continue;
        if ((qualval > 0) && (var.QUAL() < qualval)) continue;
        if (multiallelics && (!var.isMultiAllelics())) continue;
        if (multisnps && (!var.isMultiAllelicSNP())) continue;
        if (snps && (!var.isSNP())) continue;
        if (indels && (!var.isIndel())) continue;
        pos.push_back(var.POS());
        qual.push_back(var.QUAL());
        chr.push_back(var.CHROM());
        id.push_back(var.ID());
        ref.push_back(var.REF());
        alt.push_back(var.ALT());
        filter.push_back(var.FILTER());
        if (INFO) info.push_back(var.INFO());
        var.getGenotypes(gt);
        GT.push_back(gt);
    }
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

// [[Rcpp::export]]
List tableOther(std::string format, std::string vcffile, std::string region,
                std::string samples = "-", double qualval = 0, bool pass = false, bool INFO = true,
                bool snps = false, bool indels = false, bool multiallelics = false,
                bool multisnps = false) {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<int> pos;
    vector<float> qual;
    vector<std::string> chr, ref, alt, id, filter, info;
    int tagtype = vcf.header.getFormatType(format);
    if (tagtype == 1) {
        vector<vector<int>> mat;
        vector<int> vec;
        while (vcf.getNextVariant(var)) {
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.INFO());
            var.getFORMAT(format, vec);
            int nvals = vec.size() / vcf.nsamples;  // how many values per sample
            for (int i = 0; i < vcf.nsamples; i++) {
                for (int j = 0; j < nvals; j++)
                    // hit the end, set it to NA == bcf_int32_missing
                    if (vec[i * nvals + j] == bcf_int32_vector_end ||
                        vec[i * nvals + j] == bcf_int32_missing)
                        vec[i * nvals + j] = NA_INTEGER;
            }
            mat.push_back(vec);
        }
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named(format) = mat);
    } else if (tagtype == 2) {
        vector<vector<double>> mat;
        vector<float> vecf;
        vector<double> vecd;
        while (vcf.getNextVariant(var)) {
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.INFO());
            var.getFORMAT(format, vecf);
            int nvals = vecf.size() / vcf.nsamples;  // how many values per sample
            vecd.resize(vecf.size());
            for (int i = 0; i < vcf.nsamples; i++) {
                for (int j = 0; j < nvals; j++)
                    if (bcf_float_is_missing(vecf[i * nvals + j]) ||
                        bcf_float_is_vector_end(vecf[i * nvals + j]))
                        vecd[i * nvals + j] = NA_REAL;
                    else
                        vecd[i * nvals + j] = vecf[i * nvals + j];
            }
            mat.push_back(vecd);
        }
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named(format) = mat);
    } else if (tagtype == 3) {
        vector<vector<std::string>> mat;
        vector<std::string> vec;
        while (vcf.getNextVariant(var)) {
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.INFO());
            var.getFORMAT(format, vec);
            mat.push_back(vec);
        }
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named(format) = mat);
    } else {
        while (vcf.getNextVariant(var)) {
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.INFO());
        }
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named("na") = IntegerVector::create(NA_INTEGER));
    }
}
