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
    vector<int> snp_count(vcf.nsamples, 0), sv_count(vcf.nsamples, 0), indel_count(vcf.nsamples, 0),
        multi_count(vcf.nsamples, 0);
    vector<int> gt;
    while (vcf.getNextVariant(var)) {
        var.getGenotypes(gt);
        for (int i = 0; i < vcf.nsamples; i++) {
            if (!var.isGenoMissing[i]) {
                if (var.isSNP()) snp_count[i] += 1;
                if (var.isIndel()) indel_count[i] += 1;
                if (var.isSV()) sv_count[i] += 1;
                if (var.isMultiAllelic()) multi_count[i] += 1;
            }
        }
    }
    return List::create(Named("SNP") = snp_count, Named("INDEL") = indel_count,
                        Named("SV") = sv_count, Named("MultiAllelic") = multi_count);
}
