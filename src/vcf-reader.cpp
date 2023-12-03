#include <Rcpp.h>
#include "vcfpp.h"

using namespace std;

class vcfreader {
   public:
    vcfreader(std::string vcffile, std::string region = "", std::string samples = "") {
        br.open(vcffile);
        if (!samples.empty()) br.setSamples(samples);
        if (!region.empty()) br.setRegion(region);
        var.init(br.header);
    }

    ~vcfreader() {}

    bool getNextVariant() { return br.getNextVariant(var); }
    vector<int> getGenotypes() {
        var.getGenotypes(gt);
        return gt;
    }

   private:
    vcfpp::BcfReader br;
    vcfpp::BcfRecord var;
    vector<int> gt;
};

RCPP_EXPOSED_CLASS(vcfreader)
RCPP_MODULE(vcfreader) {
    using namespace Rcpp;
    class_<vcfreader>("vcfreader")
        .constructor<std::string, std::string, std::string>(
            "vcfreader given vcffile, region and samples")
        .method("getvariant", &vcfreader::getNextVariant, "get next variant record")
        .method("genotypes", &vcfreader::getGenotypes, "get genotypes");
}
