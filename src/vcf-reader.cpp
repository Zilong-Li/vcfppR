#include <Rcpp.h>
#include "vcfpp.h"

using namespace std;

//' @name vcfreader
//' @title API for reading the VCF/BCF.
//' @description Type the name of the class to see its methods
//' @field new Constructor given a vcf file \itemize{
//' \item Parameter: vcffile - The path of a vcf file
//' }
//' @field new Constructor given a vcf file and the region \itemize{
//' \item Parameter: vcffile - The path of a vcf file
//' \item Parameter: region - The region to be constrained
//' }
//' @field new Constructor given a vcf file, the region and the samples \itemize{
//' \item Parameter: vcffile - The path of a vcf file
//' \item Parameter: region - The region to be constrained
//' \item Parameter: samples - The samples to be constrained. Comma separated list of samples to include (or exclude with "^" prefix).
//' }
//' @field variant Try to get next variant record. Return false if there are no more variants or hit the end of file, otherwise return true.
//' @field chr Return the CHROM field of current variant
//' @field pos Return the POS field of current variant
//' @field id Return the CHROM field of current variant
//' @field ref Return the REF field of current variant
//' @field alt Return the ALT field of current variant
//' @field qual Return the QUAL field of current variant
//' @field filter Return the FILTER field of current variant
//' @field info Return the INFO field of current variant
//' @field infoInt Return the tag value of integer type in INFO field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in INFO}
//' @field infoFloat Return the tag value of float type in INFO field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in INFO}
//' @field infoStr Return the tag value of string type in INFO field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in INFO}
//' @field infoIntVec Return the tag value in a vector of integer type in INFO field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in INFO}
//' @field infoFloatVec Return the tag value in a vector of float type in INFO field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in INFO}
//' @field genotypes Return the genotype values in a vector of integers  \itemize{ \item Parameter: collapse - Boolean value indicates wheather to collapse the size of genotypes, eg, return diploid genotypes.}
//' @field formatInt Return the tag value of integer type for each sample in FORAMT field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in FORAMT}
//' @field formatFloat Return the tag value of float type for each sample in FORAMT field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in FORAMT}
//' @field formatStr Return the tag value of string type for each sample in FORAMT field of current variant \itemize{ \item Parameter: tag - The tag name to retrieve in FORAMT}
//' @field nsamples Return the number of samples
//' @field isSNP Test if current variant is exculsively a SNP or not
//' @field isIndel Test if current variant is exculsively a INDEL or not
//' @field isSV Test if current variant is exculsively a SV or not
//' @field isMultiAllelics Test if current variant is exculsively a Multi Allelics or not
//' @field isMultiAllelicSNP Test if current variant is exculsively a Multi Biallelics (SNPs) or not
//' @field hasSNP Test if current variant has a SNP or not
//' @field hasINDEL Test if current variant has a INDEL or not
//' @field hasINS Test if current variant has a INS or not
//' @field hasDEL Test if current variant has a DEL or not
//' @field hasMNP Test if current variant has a MNP or not
//' @field hasBND Test if current variant has a BND or not
//' @field hasOTHER Test if current variant has a OTHER or not
//' @field hasOVERLAP Test if current variant has a OVERLAP or not
class vcfreader {
   public:
    vcfreader(std::string vcffile) {
        br.open(vcffile);
        var.init(br.header);
    }

    vcfreader(std::string vcffile, std::string region) {
        br.open(vcffile);
        if (!region.empty()) br.setRegion(region);
        var.init(br.header);
    }

    vcfreader(std::string vcffile, std::string region, std::string samples) {
        br.open(vcffile);
        if (!samples.empty()) br.setSamples(samples);
        if (!region.empty()) br.setRegion(region);
        var.init(br.header);
    }

    ~vcfreader() {}

    bool variant() { return br.getNextVariant(var); }

    inline std::string chr() const { return var.CHROM(); }
    inline std::string id() const { return var.ID(); }
    inline std::string ref() const { return var.REF(); }
    inline std::string alt() const { return var.ALT(); }
    inline int pos() const { return var.POS(); }
    inline double qual() { return var.QUAL(); }
    inline std::string filter() { return var.FILTER(); }
    inline std::string info() { return var.INFO(); }

    int infoInt(std::string tag) {
        int i;
        var.getINFO(tag, i);
        return i;
    }
    double infoFloat(std::string tag) {
        float f;
        var.getINFO(tag, f);
        return (double)f;
    }
    std::string infoStr(std::string tag) {
        std::string s;
        var.getINFO(tag, s);
        return s;
    }
    vector<int> infoIntVec(std::string tag) {
        var.getINFO(tag, v_int);
        return v_int;
    }
    vector<double> infoFloatVec(std::string tag) {
        var.getINFO(tag, v_float);
        return vector<double>(v_float.begin(), v_float.end());
    }

    vector<int> genotypes(bool collapse) {
        var.getGenotypes(v_int);
        if (var.ploidy() == 2 && collapse) {
            for (size_t i = 0; i < v_int.size(); i += 2) {
                v_int[i + 1] += v_int[i];
                if (v_int[i + 1] < 0) v_int[i + 1] = NA_INTEGER;
            }
            for (size_t i = 0; 2 * i + 1 < v_int.size(); i++) {
                std::swap(v_int[i], v_int[2 * i + 1]);
            }
            v_int.resize(v_int.size() / 2);
        } else {
            for (auto& g : v_int) {
                if (g < 0) g = NA_INTEGER;
            }
        }
        return v_int;
    }

    vector<int> formatInt(std::string tag) {
        var.getFORMAT(tag, v_int);
        int nvals = v_int.size() / br.nsamples;  // how many values per sample
        for (int i = 0; i < br.nsamples; i++) {
            for (int j = 0; j < nvals; j++)
                // hit the end, set it to NA == bcf_int32_missing
                if (v_int[i * nvals + j] == bcf_int32_vector_end ||
                    v_int[i * nvals + j] == bcf_int32_missing)
                    v_int[i * nvals + j] = NA_INTEGER;
        }
        return v_int;
    }

    vector<double> formatFloat(std::string tag) {
        vector<double> vecd;
        var.getFORMAT(tag, v_float);
        int nvals = v_float.size() / br.nsamples;  // how many values per sample
        vecd.resize(v_float.size());
        for (int i = 0; i < br.nsamples; i++) {
            for (int j = 0; j < nvals; j++)
                if (bcf_float_is_missing(v_float[i * nvals + j]) ||
                    bcf_float_is_vector_end(v_float[i * nvals + j]))
                    vecd[i * nvals + j] = NA_REAL;
                else
                    vecd[i * nvals + j] = v_float[i * nvals + j];
        }
        return vecd;
    }

    vector<std::string> formatStr(std::string tag) {
        var.getFORMAT(tag, v_str);
        return v_str;
    }

    inline int nsamples() const { return br.nsamples; }
    inline bool isSNP() const { return var.isSNP(); }
    inline bool isIndel() const { return var.isIndel(); }
    inline bool isSV() const { return var.isSV(); }
    inline bool isMultiAllelics() const { return var.isMultiAllelics(); }
    inline bool isMultiAllelicSNP() const { return var.isMultiAllelicSNP(); }
    inline bool hasSNP() const { return var.hasSNP(); }
    inline bool hasINDEL() const { return var.hasINDEL(); }
    inline bool hasINS() const { return var.hasINS(); }
    inline bool hasDEL() const { return var.hasDEL(); }
    inline bool hasMNP() const { return var.hasMNP(); }
    inline bool hasBND() const { return var.hasBND(); }
    inline bool hasOTHER() const { return var.hasOTHER(); }
    inline bool hasOVERLAP() const { return var.hasOVERLAP(); }

    // WRITE
    inline void output(std::string vcffile) {
        bw.open(vcffile);
        bw.initalHeader(br.header);
    }
    inline void write() { bw.writeRecord(var); }
    inline void close() { bw.close(); }
    
    inline void setCHR(std::string s) { var.setCHR(s.c_str()); }
    inline void setID(std::string s) { var.setID(s.c_str()); }
    inline void setPOS(int pos) { var.setPOS(pos); }
    inline void setRefAlt(std::string s) { var.setRefAlt(s.c_str()); }

   private:
    vcfpp::BcfReader br;
    vcfpp::BcfRecord var;
    vcfpp::BcfWriter bw;
    vector<int> v_int;
    vector<float> v_float;
    vector<std::string> v_str;
};

RCPP_EXPOSED_CLASS(vcfreader)
RCPP_MODULE(vcfreader) {
    using namespace Rcpp;
    class_<vcfreader>("vcfreader")
        .constructor<std::string>("construct vcfreader given vcf file")
        .constructor<std::string, std::string>("construct vcfreader given vcffile and region")
        .constructor<std::string, std::string, std::string>(
            "construct vcfreader given vcf file, region and samples")
        .method("variant", &vcfreader::variant, "get next variant record")
        .method("chr", &vcfreader::chr, "get CHROM")
        .method("id", &vcfreader::id, "get ID")
        .method("pos", &vcfreader::pos, "get POS")
        .method("ref", &vcfreader::ref, "get REF")
        .method("alt", &vcfreader::alt, "get ALT")
        .method("qual", &vcfreader::qual, "get QUAL")
        .method("filter", &vcfreader::filter, "get FILTER")
        .method("info", &vcfreader::info, "get INFO as a string")
        .method("nsamples", &vcfreader::nsamples, "get the number of smaples")
        .method("genotypes", &vcfreader::genotypes, "get genotypes")
        .method("infoInt", &vcfreader::infoInt, "get tag value of int type in INFO")
        .method("infoIntVec", &vcfreader::infoIntVec,
                "get tag values in a vector of int type in INFO")
        .method("infoFloat", &vcfreader::infoFloat, "get tag value of float type in INFO")
        .method("infoFloatVec", &vcfreader::infoFloatVec,
                "get tag values in a vector of float type in INFO")
        .method("infoStr", &vcfreader::infoStr, "get tag value of string type in INFO")
        .method("formatInt", &vcfreader::formatInt, "get tag value of int type in FORMAT")
        .method("formatFloat", &vcfreader::formatFloat, "get tag value of float type in FORMAT")
        .method("formatStr", &vcfreader::formatStr, "get tag value of string type in FORMAT")
        .method("isSNP", &vcfreader::isSNP)
        .method("isIndel", &vcfreader::isIndel)
        .method("isSV", &vcfreader::isSV)
        .method("isMultiAllelics", &vcfreader::isMultiAllelics)
        .method("isMultiAllelicSNP", &vcfreader::isMultiAllelicSNP)
        .method("hasSNP", &vcfreader::hasSNP)
        .method("hasINDEL", &vcfreader::hasINDEL)
        .method("hasINS", &vcfreader::hasINS)
        .method("hasDEL", &vcfreader::hasDEL)
        .method("hasMNP", &vcfreader::hasMNP)
        .method("hasBND", &vcfreader::hasBND)
        .method("hasOTHER", &vcfreader::hasOTHER)
        .method("hasOVERLAP", &vcfreader::hasOVERLAP);
}
