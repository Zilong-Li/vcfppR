Fix the following issues.

1. Removed the redundant "in R" in the title.
2. Wrote 'C++', 'htslib' 'vcfpp.h' in single quotes
3. Updated the .Rd files by adding \value and \examples for exported methods and classes
4. Added the authors of htslib as 'cph' in ‘Authors@R’ 

Comments about documenting the C++ class:

I've documented each method of the C++ class by using a nested structure inside the \section(Fields), which is the only way I found that can document the C++ class by doxygen2 automatically.
