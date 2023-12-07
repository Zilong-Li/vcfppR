#include <Rcpp.h>
#include "vcfpp.h"

using namespace std;

//' @name vcfwriter
//' @title API for writing the VCF/BCF.
//' @description Type the name of the class to see its methods
//' @field new Constructor given a vcf file \itemize{
//' \item Parameter: vcffile - The path of a vcf file. don't start with "~"
//' \item Parameter: version - The version of VCF specification
//' }
class vcfwriter {
   public:
    vcfwriter(std::string vcffile, std::string version) : bw(vcffile, version) {}
    ~vcfwriter() {}

    inline void close() { bw.close(); }

    // WRITE ONLY
    inline void writeline(std::string line) { bw.writeLine(line); }
    inline void addLine(std::string str) { bw.header.addLine(str); }
    inline void addSample(std::string str) { bw.header.addSample(str); }
    inline void addContig(std::string str) { bw.header.addContig(str); }
    inline void addFILTER(std::string id, std::string desc) { bw.header.addFILTER(id, desc); }
    inline void addINFO(std::string id, std::string number, std::string type, std::string desc) {
        bw.header.addINFO(id, number, type, desc);
    }
    inline void addFORMAT(std::string id, std::string number, std::string type, std::string desc) {
        bw.header.addFORMAT(id, number, type, desc);
    }

   private:
    vcfpp::BcfWriter bw;
};

RCPP_EXPOSED_CLASS(vcfwriter)
RCPP_MODULE(vcfwriter) {
    using namespace Rcpp;
    class_<vcfwriter>("vcfwriter")
        .constructor<std::string, std::string>("construct vcfwriter given vcf file and the version")
        .method("close", &vcfwriter::close)
        .method("writeline", &vcfwriter::writeline)
        .method("addLine", &vcfwriter::addLine)
        .method("addSample", &vcfwriter::addSample)
        .method("addContig", &vcfwriter::addContig)
        .method("addFILTER", &vcfwriter::addFILTER)
        .method("addINFO", &vcfwriter::addINFO)
        .method("addFORMAT", &vcfwriter::addFORMAT)
        ;
}
