#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List heterozygosity(std::string vcffile, std::string region = "", std::string samples = "-",
                    bool pass = false, double qual = 0) {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<int> gt;
    vector<int> hetsum(vcf.nsamples, 0);  // store the het counts
    while (vcf.getNextVariant(var)) {
        if (pass && (var.FILTER() != "PASS")) continue;
        if ((qual > 0) && (var.QUAL() < qual)) continue;
        var.getGenotypes(gt);
        // analyze SNP variant
        if (!var.isSNP()) continue;
        assert(var.ploidy() == 2);  // make sure it is diploidy
        for (int i = 0; i < gt.size() / 2; i++)
            hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }

    return List::create(Named("samples") = vcf.header.getSamples(), Named("hets") = hetsum);
}
