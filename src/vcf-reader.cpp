#include <Rcpp.h>
#include "vcfpp.h"

using namespace std;

//' @name vcfreader
//' @title API for manipulating the VCF/BCF.
//' @description Type the name of the class to see the details and methods
//' @return A C++ class with the following fields/methods for manipulating the VCF/BCF
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
//' @field variant Try to get next variant record. return FALSE if there are no more variants or hit the end of file, otherwise TRUE.
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
//' @field nsamples Return the number of samples
//' @field samples Return a vector of samples id
//' @field header Return the raw string of the vcf header
//' @field string Return the raw string of current variant including newline
//' @field line Return the raw string of current variant without newline
//' @field output Init an output object for streaming out the variants to another vcf
//' @field write Streaming out current variant the output vcf
//' @field close Close the connection to the output vcf
//' @field setCHR Modify the CHR of current variant \itemize{ \item Parameter: s - A string for CHR}
//' @field setID Modify the ID of current variant \itemize{ \item Parameter: s - A string for ID}
//' @field setPOS Modify the POS of current variant \itemize{ \item Parameter: pos - An integer for POS}
//' @field setRefAlt Modify the REF and ALT of current variant \itemize{ \item Parameter: s - A string reperated by comma}
//' @field setInfoInt Modify the given tag of INT type in the INFO of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: v - An integer for the tag value}
//' @field setInfoFloat Modify the given tag of FLOAT type in the INFO of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: v - A double for the tag value}
//' @field setInfoStr Modify the given tag of STRING type in the INFO of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: s - A string for the tag value}
//' @field setPhasing Modify the phasing status of each sample
//' \itemize{\item Parameter: v - An integer vector with size of the number of samples. only 1s and 0s are valid.}
//' @field setGenotypes Modify the genotypes of current variant
//' \itemize{\item Parameter: v - An integer vector for genotypes. Use NA or -9 for missing value.}
//' @field setFormatInt Modify the given tag of INT type in the FORMAT of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: v - An integer for the tag value}
//' @field setFormatFloat Modify the given tag of FLOAT type in the FORMAT of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: v - A double for the tag value}
//' @field setFormatStr Modify the given tag of STRING type in the FORMAT of current variant
//' \itemize{
//' \item Parameter: tag - A string for the tag name
//' \item Parameter: s - A string for the tag value}
//' @field rmInfoTag Remove the given tag from the INFO of current variant
//' \itemize{\item Parameter: s - A string for the tag name}
//' @field rmFormatTag Remove the given tag from the FORMAT of current variant
//' \itemize{\item Parameter: s - A string for the tag name}
//' @field setVariant Modify current variant by adding a vcf line
//' \itemize{\item Parameter: s - A string for one line in the VCF}
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
//' @examples
//' vcffile <- system.file("extdata", "raw.gt.vcf.gz", package="vcfppR")
//' br <- vcfreader$new(vcffile)
//' res <- rep(0L, br$nsamples())
//' while(br$variant()) {
//'   if(br$isSNP()) {
//'   gt <- br$genotypes(TRUE) == 1
//'   gt[is.na(gt)] <- FALSE
//'   res <- res + gt
//'   }
//' }
class vcfreader {
   public:
    vcfreader(const std::string& vcffile) : fin(vcffile) {
        br.open(vcffile);
        var.initHeader(br.header);
    }

    vcfreader(const std::string& vcffile, const std::string& region) : fin(vcffile) {
        br.open(vcffile);
        if (!region.empty()) br.setRegion(region);
        var.initHeader(br.header);
    }

    vcfreader(const std::string& vcffile, const std::string& region, const std::string& samples)
        : fin(vcffile) {
        br.open(vcffile);
        if (!samples.empty()) br.setSamples(samples);
        if (!region.empty()) br.setRegion(region);
        var.initHeader(br.header);
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
    inline std::string info() { return var.allINFO(); }

    int infoInt(std::string tag) {
        int i;
        if (var.getINFO(tag, i)) {
            return i;
        } else {
            return NA_INTEGER;
        }
    }
    double infoFloat(std::string tag) {
        float f;
        if (var.getINFO(tag, f)) {
            return (double)f;
        } else {
            return NA_REAL;
        }
    }
    std::string infoStr(std::string tag) {
        std::string s{""};
        var.getINFO(tag, s);
        return s;
    }
    vector<int> infoIntVec(std::string tag) {
        if (var.getINFO(tag, v_int)) {
            return v_int;
        } else {
            return vector<int>();
        }
    }
    vector<double> infoFloatVec(std::string tag) {
        if (var.getINFO(tag, v_float)) {
            return vector<double>(v_float.begin(), v_float.end());
        } else {
            return vector<double>();
        }
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
        if (!var.getFORMAT(tag, v_int)) { return vector<int>(); }
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
        if (!var.getFORMAT(tag, v_float)) return vecd;
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
        if (var.getFORMAT(tag, v_str)) {
            return v_str;
        } else {
            return vector<std::string>();
        }
    }

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
    inline int nsamples() const { return br.nsamples; }
    inline int ploidy() {
        auto v = genotypes(false);
        return v.size() / nsamples();
    }
    inline std::string header() const { return br.header.asString(); }
    inline std::vector<std::string> samples() const { return br.header.getSamples(); }
    inline std::string string() const { return var.asString(); }
    inline std::string line() {
        std::string s = var.asString();
        s.pop_back();
        return s;
    }

    // WRITE
    inline void output(const std::string& vcffile) {
        bw.open(vcffile);
        bw.copyHeader(fin);
        var.resetHeader(bw.header);
        writable = true;
    }
    inline void write() {
        if (writable) bw.writeRecord(var);
    }
    inline void close() {
        if (writable) bw.close();
    }

    inline void setCHR(std::string s) { var.setCHR(s.c_str()); }
    inline void setID(std::string s) { var.setID(s.c_str()); }
    inline void setPOS(int pos) { var.setPOS(pos); }
    inline void setRefAlt(const std::string& s) { var.setRefAlt(s.c_str()); }
    inline bool setInfoInt(std::string tag, int v) { return var.setINFO(tag, v); }
    inline bool setInfoFloat(std::string tag, double v) { return var.setINFO(tag, v); }
    inline bool setInfoStr(std::string tag, const std::string& s) { return var.setINFO(tag, s); }
    inline bool setFormatInt(std::string tag, const vector<int>& v) {
        return var.setFORMAT(tag, v);
    }
    inline bool setFormatFloat(std::string tag, const vector<double>& v) {
        vector<float> f(v.begin(), v.end());
        return var.setFORMAT(tag, f);
    }
    inline bool setFormatStr(std::string tag, const std::string& s) {
        if (s.length() % nsamples() != 0) {
            Rcpp::Rcout << "the length of s must be divisable by nsamples()";
            return false;
        }
        return var.setFORMAT(tag, s);
    }
    inline bool setGenotypes(const vector<int>& v) {
        if ((int)v.size() != nsamples() * ploidy()) {
            Rcpp::Rcout << "nsamples: " << nsamples() << ", ploidy: " << ploidy() << "\n";
            Rcpp::Rcout << "the size of genotype vector is not equal to nsamples * ploidy\n";
            return false;
        }
        return var.setGenotypes(v);
    }
    inline void setPhasing(const vector<int>& v) {
        vector<char> c(v.begin(), v.end());
        var.setPhasing(c);
    }

    inline void rmInfoTag(std::string s) { var.removeINFO(s); }
    inline void rmFormatTag(std::string s) { var.removeFORMAT(s); }
    inline void addINFO(const std::string& id, const std::string& number, const std::string& type,
                        const std::string& desc) {
        if (writable)
            bw.header.addINFO(id, number, type, desc);
        else
            Rcpp::Rcout << "please call the `output(filename)` function first\n";
    }
    inline void addFORMAT(const std::string& id, const std::string& number, const std::string& type,
                          const std::string& desc) {
        if (writable)
            bw.header.addFORMAT(id, number, type, desc);
        else
            Rcpp::Rcout << "please call the `output(filename)` function first\n";
    }

   private:
    bool writable = false;
    const std::string fin;
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
        .method("hasOVERLAP", &vcfreader::hasOVERLAP)
        .method("nsamples", &vcfreader::nsamples)
        .method("samples", &vcfreader::samples)
        .method("header", &vcfreader::header)
        .method("string", &vcfreader::string)
        .method("line", &vcfreader::line)
        .method("output", &vcfreader::output)
        .method("write", &vcfreader::write)
        .method("close", &vcfreader::close)
        .method("setCHR", &vcfreader::setCHR)
        .method("setID", &vcfreader::setID)
        .method("setPOS", &vcfreader::setPOS)
        .method("setRefAlt", &vcfreader::setRefAlt)
        .method("setInfoInt", &vcfreader::setInfoInt)
        .method("setInfoFloat", &vcfreader::setInfoFloat)
        .method("setInfoStr", &vcfreader::setInfoStr)
        .method("setGenotypes", &vcfreader::setGenotypes)
        .method("setPhasing", &vcfreader::setPhasing)
        .method("setFormatInt", &vcfreader::setFormatInt)
        .method("setFormatFloat", &vcfreader::setFormatFloat)
        .method("setFormatStr", &vcfreader::setFormatStr)
        .method("addINFO", &vcfreader::addINFO)
        .method("addFORMAT", &vcfreader::addFORMAT)
        .method("rmInfoTag", &vcfreader::rmInfoTag)
        .method("rmFormatTag", &vcfreader::rmFormatTag);
}
