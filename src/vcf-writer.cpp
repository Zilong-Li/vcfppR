#include <Rcpp.h>
#include "vcfpp.h"

using namespace std;

//' @name vcfwriter
//' @title API for writing the VCF/BCF.
//' @description Type the name of the class to see the details and methods
//' @return A C++ class with the following fields/methods for writing the VCF/BCF
//' @field new Constructor given a vcf file \itemize{
//' \item Parameter: vcffile - The path of a vcf file. don't start with "~"
//' \item Parameter: version - The version of VCF specification
//' }
//' @field addContig Add a Contig in the header of the vcf
//' \itemize{ \item Parameter: str - A string for the CONTIG name }
//' @field addFILTER Add a FILTER in the header of the vcf
//' \itemize{
//' \item Parameter: id - A string for the FILTER name
//' \item Parameter: desc - A string for description of what it means}
//' @field addINFO Add a INFO in the header of the vcf
//' \itemize{
//' \item Parameter: id - A string for the tag name
//' \item Parameter: number - A string for the number
//' \item Parameter: type - A string for the type
//' \item Parameter: desc - A string for description of what it means}
//' @field addFORMAT Add a FORMAT in the header of the vcf
//' \itemize{
//' \item Parameter: id - A string for the tag name
//' \item Parameter: number - A string for the number
//' \item Parameter: type - A string for the type
//' \item Parameter: desc - A string for description of what it means}
//' @field addSample Add a SAMPLE in the header of the vcf
//' \itemize{ \item Parameter: str - A string for a SAMPLE name }
//' @field addLine Add a line in the header of the vcf
//' \itemize{ \item Parameter: str - A string for a line in the header of VCF }
//' @field writeline Write a variant record given a line
//' \itemize{ \item Parameter: line - A string for a line in the variant of VCF. Not ended with "newline" }
//' @field close Close and save the vcf file
//' @examples
//' outvcf <- file.path(paste0(tempfile(), ".vcf.gz"))
//' bw <- vcfwriter$new(outvcf, "VCF4.1")
//' bw$addContig("chr20")
//' bw$addFORMAT("GT", "1", "String", "Genotype");
//' bw$addSample("NA12878")
//' s1 <- "chr20\t2006060\t.\tG\tC\t100\tPASS\t.\tGT\t1|0"
//' bw$close()
class vcfwriter {
   public:
    vcfwriter(std::string vcffile, std::string version) : bw(vcffile, version) {}
    ~vcfwriter() {}

    inline void close() { bw.close(); }

    // WRITE ONLY
    inline void addContig(const std::string& str) { bw.header.addContig(str); }
    inline void addFILTER(const std::string& id, const std::string& desc) {
        bw.header.addFILTER(id, desc);
    }
    inline void addINFO(const std::string& id, const std::string& number, const std::string& type,
                        const std::string& desc) {
        bw.header.addINFO(id, number, type, desc);
    }
    inline void addFORMAT(const std::string& id, const std::string& number, const std::string& type,
                          const std::string& desc) {
        bw.header.addFORMAT(id, number, type, desc);
    }
    inline void addSample(const std::string & str) { bw.header.addSample(str); }
    inline void addLine(const std::string & str) { bw.header.addLine(str); }
    inline void writeline(const std::string & line) { bw.writeLine(line); }

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
        .method("addFORMAT", &vcfwriter::addFORMAT);
}
