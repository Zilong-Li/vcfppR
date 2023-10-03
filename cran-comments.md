## R CMD check results


❯ checking top-level files ... WARNING
  A complete check needs the 'checkbashisms' script.
  See section ‘Configure and cleanup’ in the ‘Writing R Extensions’
  manual.

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

❯ checking pragmas in C/C++ headers and code ... NOTE
  File which contains pragma(s) suppressing diagnostics:
    ‘src/htslib-1.18/htscodecs/htscodecs/fqzcomp_qual.c’

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mavx2’ ‘-mavx512f’ ‘-mpopcnt’

0 errors ✔ | 1 warning ✖ | 3 notes ✖
