#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

//' report the stats of variants
//' @param vcffile path to the VCF file with index
//' @param region  region to extract, default "" for all
//' @param samples samples to extract, default "-" for all
//' @return the counts of each type of variant
//' @export
// [[Rcpp::export]]
List summaryVariants(std::string vcffile, std::string region = "", std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    int nsamples = vcf.nsamples;
    vector<int> gt;
    vector<int> snp_count(nsamples, 0), sv_count(nsamples, 0), indel_count(nsamples, 0),
        multi_count(nsamples, 0), ins_count(nsamples, 0), del_count(nsamples, 0),
        mnp_count(nsamples, 0), multisnp_count(nsamples, 0), other_count(nsamples, 0);
    int snp{0}, indel{0}, ins{0}, del{0}, sv{0}, multiallelic{0}, mnp{0}, all{0}, multisnp{0},
        other{0};
    bool hassnp, hasmnp, hasindel, hasins, hasdel, hasother, ismulti, issv, issnpmulti;
    while (vcf.getNextVariant(var)) {
        all++;
        if ((ismulti = var.isMultiAllelics())) multiallelic += 1;
        if ((issnpmulti = var.isMultiAllelicSNP())) multisnp += 1;
        if (ismulti) continue;
        if ((hassnp = var.hasSNP())) snp += 1;
        if ((hasmnp = var.hasMNP())) mnp += 1;
        if ((hasindel = var.hasINDEL())) indel += 1;
        if ((hasins = var.hasINS())) ins += 1;
        if ((hasdel = var.hasDEL())) del += 1;
        if ((hasother = var.hasOTHER())) other += 1;
        if ((issv = var.isSV())) sv += 1;
        var.getGenotypes(gt);
        for (int i = 0; i < vcf.nsamples; i++) {
            if (!var.isGenoMissing[i]) {
                if (hassnp) snp_count[i] += 1;
                if (hasmnp) mnp_count[i] += 1;
                if (hasindel) indel_count[i] += 1;
                if (hasins) ins_count[i] += 1;
                if (hasdel) del_count[i] += 1;
                if (hasother) other_count[i] += 1;
                if (issv) sv_count[i] += 1;
            }
        }
    }

    IntegerVector stats = IntegerVector::create(
        Named("ALL", all), Named("SNP", snp), Named("MNP", mnp), Named("INDEL", indel),
        Named("INS", ins), Named("DEL", del), Named("SV", sv), Named("MultiAllelics", multiallelic),
        Named("MultiAllelicSNP", multisnp), Named("Other", other));

    return List::create(Named("summary") = stats, Named("samples") = vcf.header.getSamples(),
                        Named("SNP") = snp_count, Named("MNP") = mnp_count,
                        Named("INDEL") = indel_count, Named("INS") = ins_count,
                        Named("DEL") = del_count, Named("SV") = sv_count,
                        Named("Other") = other_count);
}
