## R CMD check results


❯ checking top-level files ... WARNING
  A complete check needs the 'checkbashisms' script.
  See section ‘Configure and cleanup’ in the ‘Writing R Extensions’
  manual.

❯ checking for GNU extensions in Makefiles ... WARNING
  Found the following file(s) containing GNU extensions:
    src/htslib/Makefile
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section ‘Writing portable
  packages’ in the ‘Writing R Extensions’ manual.

❯ checking pragmas in C/C++ headers and code ... NOTE
  File which contains pragma(s) suppressing diagnostics:
    ‘src/htslib/htscodecs/htscodecs/fqzcomp_qual.c’

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mavx2’ ‘-mavx512f’ ‘-mpopcnt’

0 errors ✔ | 2 warnings ✖ | 2 notes ✖
