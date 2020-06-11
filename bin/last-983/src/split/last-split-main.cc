// Copyright 2013, 2014 Martin C. Frith

// Parse the command line and run last-split.

#include "last-split.hh"

#include "stringify.hh"

#include <getopt.h>

#include <cctype>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>

static size_t defaultBytes(bool isSplicedAlignment) {
  size_t b = isSplicedAlignment ? 8 : 8 * 1024;
  for (int i = 0; i < 3; ++i) {
    size_t n = b * 1024;
    if (n / 1024 != b) return -1;
    b = n;
  }
  return b;
}

static char parseOutputFormat(const char *text) {
  std::string s = text;
  for (size_t i = 0; i < s.size(); ++i) {
    s[i] = std::tolower(s[i]);
  }
  if (s == "maf")  return 'm';
  if (s == "maf+") return 'M';
  return 0;
}

static void run(int argc, char* argv[]) {
  LastSplitOptions opts;

  opts.format = 'M';
  opts.direction = 1;
  opts.cis = 0.004;
  opts.trans = 1e-05;
  opts.mean = 7.0;
  opts.sdev = 1.7;
  opts.mismap = 1.0;
  opts.score = -1;
  opts.no_split = false;
  opts.bytes = 0;
  opts.verbose = false;
  opts.isSplicedAlignment = false;

  std::string version = "\
last-split "
#include "version.hh"
"\n\
License GPLv3+: GNU GPL version 3 or later.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.\n\
";

  std::string help = "\
Usage: " + std::string(argv[0]) + " [options] LAST-alignments.maf\n\
\n\
Read alignments of query sequences to a genome, and estimate the genomic\n\
source of each part of each query, allowing different parts of one query to\n\
come from different parts of the genome.\n\
\n\
Options:\n\
 -h, --help         show this help message and exit\n\
 -f, --format=FMT   output format: MAF, MAF+ (default=MAF+)\n\
 -g, --genome=NAME  lastdb genome name\n\
 -d, --direction=D  RNA direction: 0=reverse, 1=forward, 2=mixed (default="
    + cbrc::stringify(opts.direction) + ")\n\
 -c, --cis=PROB     cis-splice probability per base (default="
    + cbrc::stringify(opts.cis) + ")\n\
 -t, --trans=PROB   trans-splice probability per base (default="
    + cbrc::stringify(opts.trans) + ")\n\
 -M, --mean=MEAN    mean of ln[intron length] (default="
    + cbrc::stringify(opts.mean) + ")\n\
 -S, --sdev=SDEV    standard deviation of ln[intron length] (default="
    + cbrc::stringify(opts.sdev) + ")\n\
 -m, --mismap=PROB  maximum mismap probability (default="
    + cbrc::stringify(opts.mismap) + ")\n\
 -s, --score=INT    minimum alignment score (default=e OR e+t*ln[100])\n\
 -n, --no-split     write original, not split, alignments\n\
 -b, --bytes=B      maximum memory (default=8T for split, 8G for spliced)\n\
 -v, --verbose      be verbose\n\
 -V, --version      show version information and exit\n\
";

  const char sOpts[] = "hf:g:d:c:t:M:S:m:s:nb:vV";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "format",   required_argument, 0, 'f' },
    { "genome",   required_argument, 0, 'g' },
    { "direction",required_argument, 0, 'd' },
    { "cis",      required_argument, 0, 'c' },
    { "trans",    required_argument, 0, 't' },
    { "mean",     required_argument, 0, 'M' },
    { "sdev",     required_argument, 0, 'S' },
    { "mismap",   required_argument, 0, 'm' },
    { "score",    required_argument, 0, 's' },
    { "no-split", no_argument,       0, 'n' },
    { "bytes",    required_argument, 0, 'b' },
    { "verbose",  no_argument,       0, 'v' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return;
    case 'f':
      opts.format = parseOutputFormat(optarg);
      break;
    case 'g':
      opts.isSplicedAlignment = true;
      opts.genome = optarg;
      break;
    case 'd':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.direction, optarg);
      break;
    case 'c':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.cis, optarg);
      break;
    case 't':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.trans, optarg);
      break;
    case 'M':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.mean, optarg);
      break;
    case 'S':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.sdev, optarg);
      break;
    case 'm':
      cbrc::unstringify(opts.mismap, optarg);
      break;
    case 's':
      cbrc::unstringify(opts.score, optarg);
      break;
    case 'n':
      opts.no_split = true;
      break;
    case 'b':
      cbrc::unstringifySize(opts.bytes, optarg);
      break;
    case 'v':
      opts.verbose = true;
      break;
    case 'V':
      std::cout << version;
      return;
    case '?':
      throw std::runtime_error("");
    }
  }

  if (!opts.bytes) opts.bytes = defaultBytes(opts.isSplicedAlignment);

  opts.inputFileNames.assign(argv + optind, argv + argc);

  if (opts.inputFileNames.empty()) opts.inputFileNames.push_back("-");

  std::ios_base::sync_with_stdio(false);  // makes std::cin much faster!!!

  lastSplit(opts);
}

int main(int argc, char* argv[]) {
  try {
    run(argc, argv);
    if (!flush(std::cout)) throw std::runtime_error("write error");
    return EXIT_SUCCESS;
  } catch (const std::bad_alloc& e) {  // bad_alloc::what() may be unfriendly
    std::cerr << argv[0] << ": out of memory\n";
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    const char *s = e.what();
    if (*s) std::cerr << argv[0] << ": " << s << '\n';
    return EXIT_FAILURE;
  }
}
