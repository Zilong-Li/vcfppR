#include <Rcpp.h>
#include "vcfpp.h"

using namespace Rcpp;
using namespace std;

//' @title Calculate INFO score from GP after genotype imputation
//' @param GP vector of length a multiple of 3
//' @export
// [[Rcpp::export]]
double vcfpp_calc_info_persite(const std::vector<double>& GP){
  // if(GP.size() % 3 > 0) stop("GP must be a power of 3");
  int N = GP.size() / 3; // number of samples
  double eij = 0, fij = 0, a0, a1;
  for(int i = 0; i < N; i++) {
      a0 = GP[i * 3 + 1] + GP[i * 3 + 2] * 2;
      a1 = GP[i * 3 + 1] + GP[i * 3 + 2] * 4;
      eij += a0;
      fij += a1 - a0 * a0;
  }
  double eaf = eij / 2 / N;
  double thetaHat = std::lround(1e2 * eaf) / 1e2;
  double info = 0.0;
  if(thetaHat == 0 || thetaHat == 1){
    info = 1.0;
  } else if(info < 0) {
    info = 0.0;
  } else {
    info = 1.0 - fij / (2.0 * N * eaf * (1.0 - eaf));
    info = std::lround(1e3 * info) / 1e3;
    info = info > 1.0 ? 1.0 : info;
  }
  return info;
}

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
        if(!var.getGenotypes(gt)) continue;
        // analyze SNP variant
        if (!var.isSNP()) continue;
        assert(var.ploidy() == 2);  // make sure it is diploidy
        for (size_t i = 0; i < gt.size() / 2; i++)
            hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }

    return List::create(Named("samples") = vcf.header.getSamples(), Named("hets") = hetsum);
}

