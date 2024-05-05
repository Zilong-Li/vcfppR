#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

using UMapStringInt = std::unordered_map<std::string, int>;

UMapStringInt map_ids(const std::vector<std::string>& ids) {
    UMapStringInt ids_m;
    if (ids.size() && ids[0] != "") {
        for (auto& i : ids)
            ids_m[i] = 1;
    }
    return ids_m;
}

// [[Rcpp::export]]
List tableGT(std::string vcffile, std::string region, std::string samples, std::string format,
             const std::vector<std::string>& ids, double qualval, bool pass, bool INFO, bool snps,
             bool indels, bool multiallelics, bool multisnps, bool svs, int mac) {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<vector<int>> GT;
    vector<int> gt, pos;
    vector<float> qual;
    vector<std::string> chr, ref, alt, id, filter, info;
    UMapStringInt ids_m = map_ids(ids);
    while (vcf.getNextVariant(var)) {
        if (ids_m.size() && ids_m.count(var.ID()) == 0) continue;
        if (pass && (var.FILTER() != "PASS")) continue;
        if ((qualval > 0) && (var.QUAL() < qualval)) continue;
        if (multiallelics && (!var.isMultiAllelics())) continue;
        if (multisnps && (!var.isMultiAllelicSNP())) continue;
        if (snps && (!var.isSNP())) continue;
        if (indels && (!var.isIndel())) continue;
        if (svs && (!var.isSV())) continue;
        var.getGenotypes(gt);
        if (mac > 0) {
            int ac = 0;
            for (auto g : gt)
                ac += g;
            if (ac < mac) continue;
        }
        GT.push_back(gt);
        pos.push_back(var.POS());
        qual.push_back(var.QUAL());
        chr.push_back(var.CHROM());
        id.push_back(var.ID());
        ref.push_back(var.REF());
        alt.push_back(var.ALT());
        filter.push_back(var.FILTER());
        if (INFO) info.push_back(var.allINFO());
    }
    return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                        Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                        Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                        Named("info") = info, Named("gt") = GT);
}

// [[Rcpp::export]]
List tableFormat(std::string vcffile, std::string region, std::string samples, std::string format,
                 const std::vector<std::string>& ids, double qualval, bool pass, bool INFO,
                 bool snps, bool indels, bool multiallelics, bool multisnps, bool svs) {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<int> pos;
    vector<float> qual;
    vector<std::string> chr, ref, alt, id, filter, info;
    UMapStringInt ids_m = map_ids(ids);
    int tagtype = vcf.header.getFormatType(format);
    if (tagtype == 1) {
        vector<vector<int>> mat;
        vector<int> vec;
        while (vcf.getNextVariant(var)) {
            if (ids_m.size() && ids_m.count(var.ID()) == 0) continue;
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            if (svs && (!var.isSV())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.allINFO());
            var.getFORMAT(format, vec);
            int nvals = vec.size() / vcf.nsamples;  // how many values per sample
            for (int i = 0; i < vcf.nsamples; i++) {
                for (int j = 0; j < nvals; j++)
                    // hit the end, if true, the sample has smaller ploidy
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
            if (ids_m.size() && ids_m.count(var.ID()) == 0) continue;
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            if (svs && (!var.isSV())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.allINFO());
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
            if (ids_m.size() && ids_m.count(var.ID()) == 0) continue;
            if (pass && (var.FILTER() != "PASS")) continue;
            if ((qualval > 0) && (var.QUAL() < qualval)) continue;
            if (multiallelics && (!var.isMultiAllelics())) continue;
            if (multisnps && (!var.isMultiAllelicSNP())) continue;
            if (snps && (!var.isSNP())) continue;
            if (indels && (!var.isIndel())) continue;
            if (svs && (!var.isSV())) continue;
            pos.push_back(var.POS());
            qual.push_back(var.QUAL());
            chr.push_back(var.CHROM());
            id.push_back(var.ID());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
            filter.push_back(var.FILTER());
            if (INFO) info.push_back(var.allINFO());
            var.getFORMAT(format, vec);
            mat.push_back(vec);
        }
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named(format) = mat);
    } else {
        while (vcf.getNextVariant(var)) {
            if (ids_m.size() && ids_m.count(var.ID()) == 0) continue;
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
            if (INFO) info.push_back(var.allINFO());
        }
        Rcout << "there is no " << format << " in the FORMAT\n";
        return List::create(Named("samples") = vcf.header.getSamples(), Named("chr") = chr,
                            Named("pos") = pos, Named("id") = id, Named("ref") = ref,
                            Named("alt") = alt, Named("qual") = qual, Named("filter") = filter,
                            Named("info") = info, Named("na") = IntegerVector::create(NA_INTEGER));
    }
}
