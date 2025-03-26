/*  htsfile.c -- file identifier and minimal viewer.

    Copyright (C) 2014-2019 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif

enum { identify, view_headers, view_all, copy } mode = identify;
int show_headers = 1;
int verbose = 0;
int status = EXIT_SUCCESS;  /* Exit status from main */

void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) error(const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fprintf(stderr, "htsfile: ");
    vfprintf(stderr, format, args);
    if (err) fprintf(stderr, ": %s\n", strerror(err));
    else fprintf(stderr, "\n");
    fflush(stderr);
    va_end(args);
    status = EXIT_FAILURE;
}

static void view_sam(samFile *in, const char *filename)
{
    bam1_t *b = NULL;
    sam_hdr_t *hdr = NULL;
    samFile *out = NULL;

    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        errno = 0; error("reading headers from \"%s\" failed", filename);
        goto clean;
    }

    out = hts_open("-", "w");
    if (out == NULL) { error("reopening standard output failed"); goto clean; }

    if (show_headers) {
        if (sam_hdr_write(out, hdr) != 0) {
            error("writing headers to standard output failed");
            goto clean;
        }
    }

    if (mode == view_all) {
        int ret;

        b = bam_init1();
        if (b == NULL) { error("can't create record"); goto clean; }

        while ((ret = sam_read1(in, hdr, b)) >= 0) {
            if (sam_write1(out, hdr, b) < 0) {
                error("writing to standard output failed");
                goto clean;
            }
        }

        if (ret < -1) { error("reading \"%s\" failed", filename); goto clean; }
    }

 clean:
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    if (out) hts_close(out);
}

static void view_vcf(vcfFile *in, const char *filename)
{
    bcf1_t *rec = NULL;
    bcf_hdr_t *hdr = NULL;
    vcfFile *out = NULL;

    hdr = bcf_hdr_read(in);
    if (hdr == NULL) {
        errno = 0; error("reading headers from \"%s\" failed", filename);
        goto clean;
    }

    out = hts_open("-", "w");
    if (out == NULL) { error("reopening standard output failed"); goto clean; }

    if (show_headers) {
        if (bcf_hdr_write(out, hdr) != 0) {
            error("writing headers to standard output failed");
            goto clean;
        }
    }

    if (mode == view_all) {
        int ret;

        rec = bcf_init();
        if (rec == NULL) { error("can't create record"); goto clean; }

        while ((ret = bcf_read(in, hdr, rec)) >= 0) {
            if (bcf_write(out, hdr, rec) < 0) {
                error("writing to standard output failed");
                goto clean;
            }
        }

        if (ret < -1) { error("reading \"%s\" failed", filename); goto clean; }
    }

 clean:
    if (hdr) bcf_hdr_destroy(hdr);
    if (rec) bcf_destroy(rec);
    if (out) hts_close(out);
}

static void view_raw(hFILE *fp, const char *filename)
{
    int c, prev;
    for (prev = '\n'; (c = hgetc(fp)) != EOF; prev = c)
        if (isprint(c) || c == '\n' || c == '\t') putchar(c);
        else if (c == '\r') fputs("\\r", stdout);
        else if (c == '\0') fputs("\\0", stdout);
        else printf("\\x%02x", c);

    if (prev != '\n') putchar('\n');

    if (herrno(fp)) {
        errno = herrno(fp);
        error("reading \"%s\" failed", filename);
    }
}

static void copy_raw(const char *srcfilename, const char *destfilename)
{
    hFILE *src = hopen(srcfilename, "r");
    if (src == NULL) {
        error("can't open \"%s\"", srcfilename);
        return;
    }

    size_t bufsize = 1048576;
    char *buffer = malloc(bufsize);
    if (buffer == NULL) {
        error("can't allocate copy buffer");
        hclose_abruptly(src);
        return;
    }

    hFILE *dest = hopen(destfilename, "w");
    if (dest == NULL) {
        error("can't create \"%s\"", destfilename);
        hclose_abruptly(src);
        free(buffer);
        return;
    }

    ssize_t n;
    while ((n = hread(src, buffer, bufsize)) > 0)
        if (hwrite(dest, buffer, n) != n) {
            error("writing to \"%s\" failed", destfilename);
            hclose_abruptly(dest);
            dest = NULL;
            break;
        }

    if (n < 0) {
        error("reading from \"%s\" failed", srcfilename);
        hclose_abruptly(src);
        src = NULL;
    }

    if (dest && hclose(dest) < 0) error("closing \"%s\" failed", destfilename);
    if (src && hclose(src) < 0)   error("closing \"%s\" failed", srcfilename);
    free(buffer);
}


int main(int argc, char **argv)
{
    return 0;
}
