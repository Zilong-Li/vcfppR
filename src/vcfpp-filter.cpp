#include "vcfpp.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//' get DS (dosage) of each sample and skip variants given a INFO tag of Float type and a cutoff
//' @param cutoff  skip variants with tag value smaller than this cutoff
//' @param tag     tag exists in INFO column
//' @param vcffile path to the VCF file with index
//' @param region  region to extract
//' @param samples samples to extract
//' @return A list of genotype dosages for each sample along with variant information in VCF
//' @export
// [[Rcpp::export]]
List getDSskipInfoTagFloat(double cutoff, std::string tag, std::string vcffile, std::string region = "",
                            std::string samples = "-") {
    vcfpp::BcfReader vcf(vcffile, region, samples);
    vcfpp::BcfRecord var(vcf.header);
    vector<float> ds;
    vector<vector<float>> DS;
    vector<string> chr, ref, alt;
    vector<int> pos;
    float info;
    while (vcf.getNextVariant(var)) {
        var.getINFO(tag, info);
        if(info < cutoff) continue;
        var.getFORMAT("DS", ds);
        DS.push_back(ds);
        chr.push_back(var.CHROM());
        pos.push_back(var.POS());
        ref.push_back(var.REF());
        alt.push_back(var.ALT());
    }
    return List::create(Named("samples") = vcf.header.getSamples(),
                        Named("chr") = chr, Named("pos") = pos,
                        Named("ref") = ref, Named("alt") = alt,
                        Named("ds") = DS);
}
