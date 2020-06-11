// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "LastdbArguments.hh"
#include "stringify.hh"
#include "getoptUtil.hh"
#include <iostream>
#include <stdexcept>
#include <cstdlib>  // EXIT_SUCCESS
#include <cstring>  // strtok

#define ERR(x) throw std::runtime_error(x)

static void badopt( char opt, const char* arg ){
  ERR( std::string("bad option value: -") + opt + ' ' + arg );
}

static int myGetopt( int argc, char** argv, const char* optstring ){
  if( optind < argc ){
    std::string nextarg = argv[optind];
    if( nextarg == "--help"    ) return 'h';
    if( nextarg == "--version" ) return 'V';
  }
  return getopt( argc, argv, optstring );
}

using namespace cbrc;

LastdbArguments::LastdbArguments() :
  isProtein(false),
  isKeepLowercase(true),
  tantanSetting(0),
  isCaseSensitive(false),
  seedPatterns(0),
  strand(1),
  volumeSize(-1),
  indexStep(1),
  minimizerWindow(1),
  numOfThreads(1),
  subsetSeedFile(""),
  userAlphabet(""),
  minSeedLimit(0),
  bucketDepth(-1),  // means: use the default (adapts to the data)
  childTableType(0),
  isCountsOnly(false),
  verbosity(0),
  inputFormat(sequenceFormat::fasta){}

void LastdbArguments::fromArgs( int argc, char** argv, bool isOptionsOnly ){
  programName = argv[0];
  std::string usage = "\
Usage: " + std::string(programName) +
    " [options] output-name fasta-sequence-file(s)\n\
Prepare sequences for subsequent alignment with lastal.\n\
\n\
Main Options:\n\
-h, --help: show all options and their default settings, and exit\n\
-p: interpret the sequences as proteins\n\
-R: repeat-marking options (default="
    + stringify(isKeepLowercase) + stringify(tantanSetting) + ")\n\
-c: soft-mask lowercase letters (in reference *and* query sequences)\n\
-u: seeding scheme (default: YASS for DNA, else exact-match seeds)";

  std::string help = usage + "\n\
\n\
Advanced Options (default settings):\n\
-w: use initial matches starting at every w-th position in each sequence ("
    + stringify(indexStep) + ")\n\
-W: use \"minimum\" positions in sliding windows of W consecutive positions ("
    + stringify(minimizerWindow) + ")\n\
-S: strand: 0=reverse, 1=forward, 2=both ("
    + stringify(strand) + ")\n\
-s: volume size (unlimited)\n\
-Q: input format: 0=fasta or fastq-ignore,\n\
                  1=fastq-sanger, 2=fastq-solexa, 3=fastq-illumina (fasta)\n\
-P: number of parallel threads ("
    + stringify(numOfThreads) + ")\n\
-m: seed pattern\n\
-a: user-defined alphabet\n\
-i: minimum limit on initial matches per query position ("
    + stringify(minSeedLimit) + ")\n\
-b: bucket depth\n\
-C: child table type: 0=none, 1=byte-size, 2=short-size, 3=full ("
    + stringify(childTableType) + ")\n\
-x: just count sequences and letters\n\
-v: be verbose: write messages about what lastdb is doing\n\
-V, --version: show version information, and exit\n\
\n\
Report bugs to: last-align (ATmark) googlegroups (dot) com\n\
LAST home page: http://last.cbrc.jp/\n\
";

  static const char sOpts[] = "hVpR:cm:S:s:w:W:P:u:a:i:b:C:xvQ:";

  int c;
  while ((c = myGetopt(argc, argv, sOpts)) != -1) {
    switch(c){
    case 'h':
      std::cout << help;
      throw EXIT_SUCCESS;
    case 'V':
      std::cout << "lastdb "
#include "version.hh"
	"\n";
      throw EXIT_SUCCESS;
    case 'p':
      isProtein = true;
      break;
    case 'R':
      if( optarg[0] < '0' || optarg[0] > '1' ) badopt( c, optarg );
      if( optarg[1] < '0' || optarg[1] > '2' ) badopt( c, optarg );
      if( optarg[2] ) badopt( c, optarg );
      isKeepLowercase = optarg[0] - '0';
      tantanSetting = optarg[1] - '0';
      break;
    case 'c':
      isCaseSensitive = true;
      break;
    case 'm':
      seedPatterns.push_back(optarg);
      break;
    case 'S':
      unstringify( strand, optarg );
      if( strand < 0 || strand > 2 ) badopt( c, optarg );
      break;
    case 's':
      unstringifySize( volumeSize, optarg );
      break;
    case 'w':
      unstringify( indexStep, optarg );
      if( indexStep < 1 ) badopt( c, optarg );
      break;
    case 'W':
      unstringify( minimizerWindow, optarg );
      if( minimizerWindow < 1 ) badopt( c, optarg );
      break;
    case 'P':
      unstringify( numOfThreads, optarg );
      break;
    case 'u':
      subsetSeedFile = optarg;
      break;
    case 'a':
      userAlphabet = optarg;
      break;
    case 'i':
      unstringify( minSeedLimit, optarg );
      break;
    case 'b':
      unstringify( bucketDepth, optarg );
      break;
    case 'C':
      unstringify( childTableType, optarg );
      if( childTableType < 0 || childTableType > 3 ) badopt( c, optarg );
      break;
    case 'x':
      isCountsOnly = true;
      break;
    case 'v':
      ++verbosity;
      break;
    case 'Q':
      unstringify( inputFormat, optarg );
      if( inputFormat >= sequenceFormat::prb ) badopt( c, optarg );
      break;
    case '?':
      ERR( "bad option" );
    }
  }

  if( !isOptionsOnly ){
    if( optind >= argc )
      ERR( "please give me an output name and sequence file(s)\n\n" + usage );
    lastdbName = argv[optind++];
    inputStart = optind;
  }

  resetGetopt();
}

void LastdbArguments::fromLine( const std::string& line ){
  const char* delimiters = " \t";
  const char* s = line.c_str();
  std::vector<char> args( s, s + line.size() + 1 );
  std::vector<char*> argv;
  char* i = std::strtok( &args[0], delimiters );
  argv.push_back(i);
  while( i ){
    i = std::strtok( 0, delimiters );
    argv.push_back(i);
  }
  fromArgs( argv.size() - 1, &argv[0], true );
}

void LastdbArguments::fromString( const std::string& s ){
  std::string trigger = "#lastdb";
  std::istringstream iss(s);
  std::string line;
  while( getline( iss, line ) )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      fromLine( line );
}
