// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "LastalArguments.hh"
#include "stringify.hh"
#include "getoptUtil.hh"
#include <algorithm>  // max
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cmath>  // log
#include <cstring>  // strtok
#include <cstdlib>  // EXIT_SUCCESS

#define ERR(x) throw std::runtime_error(x)

static void badopt( char opt, const char* arg ){
  ERR( std::string("bad option value: -") + opt + ' ' + arg );
}

static void parseIntList(const char *in, std::vector<int> &out) {
  std::istringstream s(in);
  out.clear();
  int x;
  while (s >> x) {
    out.push_back(x);
    s.ignore();
  }
  if (!s.eof()) {
    ERR(std::string("bad option value: ") + in);
  }
}

static void writeIntList(std::ostream &s, const std::vector<int> &v) {
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) s << ',';
    s << v[i];
  }
}

static int myGetopt( int argc, char** argv, const char* optstring ){
  if( optind < argc ){
    std::string nextarg = argv[optind];
    if( nextarg == "--help"    ) return 'h';
    if( nextarg == "--version" ) return 'V';
  }
  return getopt( argc, argv, optstring );
}

static char parseOutputFormat( const char* text ){
  std::string s = text;
  for( size_t i = 0; i < s.size(); ++i ){
    s[i] = std::tolower( s[i] );
  }
  if( s == "tab" || s == "0" ) return 't';
  if( s == "maf" || s == "1" ) return 'm';
  if( s == "blasttab" )        return 'b';
  if( s == "blasttab+" )       return 'B';
  return 0;
}

static int parseScoreAndSuffix(const std::string &s, char &suffix) {
  char c = 0;
  int x;
  std::istringstream i(s);
  if (!(i >> x) || x < 0 || (i >> c && i >> c)) {
    ERR("bad value: " + s);
  }
  suffix = c;
  return x;
}

namespace cbrc{

LastalArguments::LastalArguments() :
  outputFormat('m'),
  outputType(3),
  strand(-1),  // depends on the alphabet
  isQueryStrandMatrix(false),
  isGreedy(false),
  globality(0),
  isKeepLowercase(true),  // depends on the option used with lastdb
  tantanSetting(-1),  // depends on the option used with lastdb
  maskLowercase(-1),  // depends on the lowercase option used with lastdb
  maxEvalue(-1),
  queryLettersPerRandomAlignment(1e6),
  minScoreGapped(-1),  // depends on the alphabet
  minScoreGapless(-1),  // depends on minScoreGapped and the outputType
  matchScore(-1),  // depends on the alphabet
  mismatchCost(-1),  // depends on the alphabet
  gapPairCost(-1),  // this means: OFF
  frameshiftCost(-1),  // this means: ordinary, non-translated alignment
  matrixFile(""),
  ambiguousLetterOpt(0),
  maxDropGapped(100),  // depends on maxDropFinal
  maxDropGappedSuffix('%'),
  maxDropGapless(-1),  // depends on the score matrix
  maxDropFinal(100),  // depends on minScoreGapped
  maxDropFinalSuffix('%'),
  inputFormat(sequenceFormat::fasta),
  minHitDepth(1),
  maxHitDepth(-1),
  oneHitMultiplicity(10),
  maxGaplessAlignmentsPerQueryPosition(0),  // depends on oneHitMultiplicity
  maxAlignmentsPerQueryStrand(-1),
  cullingLimitForGaplessAlignments(0),
  cullingLimitForFinalAlignments(0),
  queryStep(1),
  minimizerWindow(0),  // depends on the reference's minimizer window
  batchSize(0),  // depends on the outputType, and voluming
  numOfThreads(1),
  maxRepeatDistance(1000),  // sufficiently conservative?
  temperature(-1),  // depends on the score matrix
  gamma(1),
  geneticCodeFile(""),
  verbosity(0){}

void LastalArguments::fromArgs( int argc, char** argv, bool optionsOnly ){
  programName = argv[0];
  std::string usage = "Usage: " + std::string(programName) +
    " [options] lastdb-name fasta-sequence-file(s)\n\
Find and align similar sequences.\n\
\n\
Cosmetic options:\n\
-h, --help: show all options and their default settings, and exit\n\
-V, --version: show version information, and exit\n\
-v: be verbose: write messages about what lastal is doing\n\
-f: output format: TAB, MAF, BlastTab, BlastTab+ (default=MAF)";

  std::string help = usage + "\n\
\n\
E-value options (default settings):\n\
-D: query letters per random alignment ("
    + stringify(queryLettersPerRandomAlignment) + ")\n\
-E: maximum expected alignments per square giga (1e+18/D/refSize/numOfStrands)\n\
\n\
Score options (default settings):\n\
-r: match score   (2 if -M, else  6 if 0<Q<5, else 1 if DNA)\n\
-q: mismatch cost (3 if -M, else 18 if 0<Q<5, else 1 if DNA)\n\
-p: match/mismatch score matrix (protein-protein: BL62, DNA-protein: BL80)\n\
-X: N/X is ambiguous in: 0=neither sequence, 1=reference, 2=query, 3=both ("
    + stringify(ambiguousLetterOpt) + ")\n\
-a: gap existence cost (DNA: 7, protein: 11, 0<Q<5: 21)\n\
-b: gap extension cost (DNA: 1, protein:  2, 0<Q<5:  9)\n\
-A: insertion existence cost (a)\n\
-B: insertion extension cost (b)\n\
-c: unaligned residue pair cost (off)\n\
-F: frameshift cost (off)\n\
-x: maximum score drop for preliminary gapped alignments (z)\n\
-y: maximum score drop for gapless alignments (min[t*10, x])\n\
-z: maximum score drop for final gapped alignments (e-1)\n\
-d: minimum score for gapless alignments (min[e, t*ln(1000*refSize/n)])\n\
-e: minimum score for gapped alignments\n\
\n\
Initial-match options (default settings):\n\
-m: maximum initial matches per query position ("
    + stringify(oneHitMultiplicity) + ")\n\
-l: minimum length for initial matches ("
    + stringify(minHitDepth) + ")\n\
-L: maximum length for initial matches (infinity)\n\
-k: use initial matches starting at every k-th position in each query ("
    + stringify(queryStep) + ")\n\
-W: use \"minimum\" positions in sliding windows of W consecutive positions\n\
\n\
Miscellaneous options (default settings):\n\
-s: strand: 0=reverse, 1=forward, 2=both (2 for DNA, 1 for protein)\n\
-S: score matrix applies to forward strand of: 0=reference, 1=query ("
    + stringify(isQueryStrandMatrix) + ")\n\
-K: omit alignments whose query range lies in >= K others with > score (off)\n\
-C: omit gapless alignments in >= C others with > score-per-length (off)\n\
-P: number of parallel threads ("
    + stringify(numOfThreads) + ")\n\
-i: query batch size (8 KiB, unless there is > 1 thread or lastdb volume)\n\
-M: find minimum-difference alignments (faster but cruder)\n\
-T: type of alignment: 0=local, 1=overlap ("
    + stringify(globality) + ")\n\
-n: maximum gapless alignments per query position (infinity if m=0, else m)\n\
-N: stop after the first N alignments per query strand\n\
-R: repeat-marking options (the same as was used for lastdb)\n\
-u: mask lowercase during extensions: 0=never, 1=gapless,\n\
    2=gapless+postmask, 3=always (2 if lastdb -c and Q<5, else 0)\n\
-w: suppress repeats inside exact matches, offset by <= this distance ("
    + stringify(maxRepeatDistance) + ")\n\
-G: genetic code file\n\
-t: 'temperature' for calculating probabilities (1/lambda)\n\
-g: 'gamma' parameter for gamma-centroid and LAMA ("
    + stringify(gamma) + ")\n\
-j: output type: 0=match counts, 1=gapless, 2=redundant gapped, 3=gapped,\n\
                 4=column ambiguity estimates, 5=gamma-centroid, 6=LAMA,\n\
                 7=expected counts ("
    + stringify(outputType) + ")\n\
-Q: input format: 0=fasta or fastq-ignore, 1=fastq-sanger, 2=fastq-solexa,\n\
                  3=fastq-illumina, 4=prb, 5=PSSM (fasta)\n\
\n\
Report bugs to: last-align (ATmark) googlegroups (dot) com\n\
LAST home page: http://last.cbrc.jp/\n\
";

  int c;
  const char optionString[] =
    "hVvf:"
    "r:q:p:X:a:b:A:B:c:F:x:y:z:d:e:"
    "D:E:"
    "s:S:MT:m:l:L:n:N:C:K:k:W:i:P:R:u:w:t:g:G:j:Q:";
  while( (c = myGetopt(argc, argv, optionString)) != -1 ){
    switch(c){
    case 'h':
      std::cout << help;
      throw EXIT_SUCCESS;
    case 'V':
      std::cout << "lastal "
#include "version.hh"
	"\n";
      throw EXIT_SUCCESS;
    case 'v':
      ++verbosity;
      break;
    case 'f':
      outputFormat = parseOutputFormat( optarg );
      if( !outputFormat ) badopt( c, optarg );
      break;

    case 'r':
      unstringify( matchScore, optarg );
      if( matchScore <= 0 ) badopt( c, optarg );
      break;
    case 'q':
      unstringify( mismatchCost, optarg );
      if( mismatchCost < 0 ) badopt( c, optarg );  // allow 0 for Fujibuchi-san
      break;
    case 'p':
      matrixFile = optarg;
      break;
    case 'X':
      unstringify(ambiguousLetterOpt, optarg);
      if (ambiguousLetterOpt < 0 || ambiguousLetterOpt > 3) badopt(c, optarg);
      break;
    case 'a':
      parseIntList(optarg, delOpenCosts);
      break;
    case 'b':
      parseIntList(optarg, delGrowCosts);
      break;
    case 'A':
      parseIntList(optarg, insOpenCosts);
      break;
    case 'B':
      parseIntList(optarg, insGrowCosts);
      break;
    case 'c':
      unstringify( gapPairCost, optarg );
      if( gapPairCost <= 0 ) badopt( c, optarg );
      break;
    case 'F':
      unstringify( frameshiftCost, optarg );
      if( frameshiftCost < 0 ) badopt( c, optarg );
      break;
    case 'x':
      maxDropGapped = parseScoreAndSuffix(optarg, maxDropGappedSuffix);
      break;
    case 'y':
      unstringify( maxDropGapless, optarg );
      if( maxDropGapless < 0 ) badopt( c, optarg );
      break;
    case 'z':
      maxDropFinal = parseScoreAndSuffix(optarg, maxDropFinalSuffix);
      break;
    case 'd':
      unstringify( minScoreGapless, optarg );
      if( minScoreGapless < 0 ) badopt( c, optarg );
      break;
    case 'e':
      unstringify( minScoreGapped, optarg );
      if( minScoreGapped < 0 ) badopt( c, optarg );
      break;

    case 'D':
      unstringify( queryLettersPerRandomAlignment, optarg );
      if( queryLettersPerRandomAlignment <= 0 ) badopt( c, optarg );
      break;
    case 'E':
      unstringify( maxEvalue, optarg );
      if( maxEvalue <= 0 ) badopt( c, optarg );
      break;

    case 'm':
      unstringify( oneHitMultiplicity, optarg );
      break;
    case 'l':
      unstringify( minHitDepth, optarg );
      break;
    case 'L':
      unstringify( maxHitDepth, optarg );
      break;
    case 'k':
      unstringify( queryStep, optarg );
      if( queryStep <= 0 ) badopt( c, optarg );
      break;
    case 'W':
      unstringify( minimizerWindow, optarg );  // allow 0, meaning "default"
      break;

    case 's':
      unstringify( strand, optarg );
      if( strand < 0 || strand > 2 ) badopt( c, optarg );
      break;
    case 'S':
      unstringify( isQueryStrandMatrix, optarg );
      break;
    case 'M':
      isGreedy = true;
      break;
    case 'T':
      unstringify( globality, optarg );
      if( globality < 0 || globality > 1 ) badopt( c, optarg );
      break;
    case 'n':
      unstringify( maxGaplessAlignmentsPerQueryPosition, optarg );
      if( maxGaplessAlignmentsPerQueryPosition <= 0 ) badopt( c, optarg );
      break;
    case 'N':
      unstringify( maxAlignmentsPerQueryStrand, optarg );
      break;
    case 'C':
      unstringify( cullingLimitForGaplessAlignments, optarg );
      break;
    case 'K':
      unstringify( cullingLimitForFinalAlignments, optarg );
      break;
    case 'i':
      unstringifySize( batchSize, optarg );
      if( batchSize <= 0 ) badopt( c, optarg );  // 0 means "not specified"
      break;
    case 'P':
      unstringify( numOfThreads, optarg );
      break;
    case 'R':
      if( optarg[0] < '0' || optarg[0] > '1' ) badopt( c, optarg );
      if( optarg[1] < '0' || optarg[1] > '2' ) badopt( c, optarg );
      if( optarg[2] ) badopt( c, optarg );
      isKeepLowercase = optarg[0] - '0';
      tantanSetting = optarg[1] - '0';
      break;
    case 'u':
      unstringify( maskLowercase, optarg );
      if( maskLowercase < 0 || maskLowercase > 3 ) badopt( c, optarg );
      break;
    case 'w':
      unstringify( maxRepeatDistance, optarg );
      break;
    case 't':
      unstringify( temperature, optarg );
      if( temperature <= 0 ) badopt( c, optarg );
      break;
    case 'g':
      unstringify( gamma, optarg );
      if( gamma <= 0 ) badopt( c, optarg );
      break;
    case 'G':
      geneticCodeFile = optarg;
      break;
    case 'j':
      unstringify( outputType, optarg );
      if( outputType < 0 || outputType > 7 ) badopt( c, optarg );
      break;
    case 'Q':
      unstringify( inputFormat, optarg );
      break;

    case '?':
      ERR( "bad option" );
    }
  }

  if( maskLowercase == 1 && inputFormat == 5 )
    ERR( "can't combine option -u 1 with option -Q 5" );

  if( maskLowercase == 2 && inputFormat == 5 )
    ERR( "can't combine option -u 2 with option -Q 5" );

  if( isTranslated() && inputFormat == 5 )
    ERR( "can't combine option -F with option -Q 5" );

  if( isFrameshift() && outputType > 3 )
    ERR( "can't combine option -F > 0 with option -j > 3" );

  if( isFrameshift() && globality == 1 )
    ERR( "can't combine option -F > 0 with option -T 1" );

  if( isTranslated() && isQueryStrandMatrix )
    ERR( "can't combine option -F with option -S 1" );

  if( isGreedy && outputType > 3 )
    ERR( "can't combine option -M with option -j > 3" );

  if( isGreedy && globality == 1 )
    ERR( "can't combine option -M with option -T 1" );

  if( isGreedy && maskLowercase == 3 )
    ERR( "can't combine option -M with option -u 3" );

  if( !optionsOnly ){
    if( optind >= argc )
      ERR( "please give me a database name and sequence file(s)\n\n" + usage );
    lastdbName = argv[optind++];
    inputStart = optind;
  }

  resetGetopt();
}

void LastalArguments::fromLine( const std::string& line ){
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

void LastalArguments::fromString( const std::string& s ){
  std::string trigger = "#last";
  std::istringstream iss(s);
  std::string line;
  while( std::getline( iss, line ) )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      fromLine( line );
}

const char* LastalArguments::matrixName( bool isProtein ) const{
  if( matrixFile.empty() && matchScore < 0 && mismatchCost < 0 && isProtein )
    return isTranslated() ? "BL80" : "BL62";
  return matrixFile.c_str();
}

void LastalArguments::setDefaultsFromAlphabet( bool isDna, bool isProtein,
					       double numLettersInReference,
					       bool isKeepRefLowercase,
					       int refTantanSetting,
                                               bool isCaseSensitiveSeeds,
					       bool isVolumes,
					       size_t refMinimizerWindow,
					       unsigned realNumOfThreads ){
  if( strand < 0 ) strand = (isDna || isTranslated()) ? 2 : 1;

  if( isGreedy ){
    if( matchScore     < 0 ) matchScore     =   2;
    if( mismatchCost   < 0 ) mismatchCost   =   3;
    int gapGrowCost = mismatchCost + matchScore / 2;
    delOpenCosts.assign(1, 0);
    delGrowCosts.assign(1, gapGrowCost);
    insOpenCosts.assign(1, 0);
    insGrowCosts.assign(1, gapGrowCost);
    if( frameshiftCost > 0 ) frameshiftCost =   0;
    if( matchScore % 2 )
      ERR( "with option -M, the match score (option -r) must be even" );
  }
  else if( isProtein ){
    // default match & mismatch scores: Blosum62 matrix
    if( matchScore < 0 && mismatchCost >= 0 ) matchScore   = 1;  // idiot-proof
    if( mismatchCost < 0 && matchScore >= 0 ) mismatchCost = 1;  // idiot-proof
    if (delOpenCosts.empty()) delOpenCosts.assign(1, 11);
    if (delGrowCosts.empty()) delGrowCosts.assign(1,  2);
  }
  else if( !isUseQuality( inputFormat ) ){
    if( matchScore     < 0 ) matchScore     =   1;
    if( mismatchCost   < 0 ) mismatchCost   =   1;
    if (delOpenCosts.empty()) delOpenCosts.assign(1, 7);
    if (delGrowCosts.empty()) delGrowCosts.assign(1, 1);
  }
  else{  // sequence quality scores will be used:
    if( matchScore     < 0 ) matchScore     =   6;
    if( mismatchCost   < 0 ) mismatchCost   =  18;
    if (delOpenCosts.empty()) delOpenCosts.assign(1, 21);
    if (delGrowCosts.empty()) delGrowCosts.assign(1,  9);
    // With this scoring scheme for DNA, gapless lambda ~= ln(10)/10,
    // so these scores should be comparable to PHRED scores.
    // Furthermore, since mismatchCost/matchScore = 3, the target
    // distribution of paired bases ~= 99% identity.  Because the
    // quality scores are unlikely to be perfect, it may be best to
    // use a lower target %identity than we otherwise would.
  }

  if( maxEvalue < 0 ){
    double r = std::max( numLettersInReference, 1.0 );  // avoid divide-by-0
    maxEvalue = 1e18 / (numOfStrands() * r * queryLettersPerRandomAlignment);
  }

  if (insOpenCosts.empty()) insOpenCosts = delOpenCosts;
  if (insGrowCosts.empty()) insGrowCosts = delGrowCosts;

  if( tantanSetting < 0 ){
    isKeepLowercase = isKeepRefLowercase;
    tantanSetting = refTantanSetting;
  }

  if( maskLowercase < 0 ){
    if( isCaseSensitiveSeeds && inputFormat != sequenceFormat::pssm )
      maskLowercase = 2;
    else
      maskLowercase = 0;
  }

  if( batchSize == 0 ){
    // With voluming, we want the batches to be as large as will
    // comfortably fit into memory, because each volume gets read from
    // disk once per batch.  With multi-threads, we want large batches
    // so that long query sequences can be processed in parallel.
    if( !isVolumes && realNumOfThreads == 1 )
      batchSize = 0x2000;  // 8 Kbytes (?)
    else if( inputFormat == sequenceFormat::pssm )
      batchSize = 0x100000;   // 1 Mbyte
    else if( outputType == 0 )
      batchSize = 0x1000000;  // 16 Mbytes
    else if( !isVolumes )
      batchSize = 0x800000;   // 8 Mbytes
    else if( inputFormat == sequenceFormat::prb )
      batchSize = 0x2000000;  // 32 Mbytes (?)
    else
      batchSize = 0x4000000;  // 64 Mbytes (?)
    // 128 Mbytes seemed to sometimes use excessive memory to store
    // the alignments.  I suspect 64 Mbytes may still be too much
    // sometimes.
    if( verbosity )
      std::cerr << programName << ": batch size=" << batchSize << '\n';
  }

  if( maxGaplessAlignmentsPerQueryPosition == 0 )
    maxGaplessAlignmentsPerQueryPosition =
      (oneHitMultiplicity > 0) ? oneHitMultiplicity : -1;

  if( minimizerWindow == 0 ) minimizerWindow = refMinimizerWindow;

  if (delOpenCosts.size() != delGrowCosts.size() ||
      insOpenCosts.size() != insGrowCosts.size()) {
    ERR("bad gap costs");
  }

  if (delOpenCosts.size() > 1 || insOpenCosts.size() > 1) {
    ERR("piecewise linear gap costs not implemented");
  }

  if (isFrameshift()) {
    if (frameshiftCost < delGrowCosts[0])
      ERR("the frameshift cost must not be less than the gap extension cost");

    if (insOpenCosts[0] != delOpenCosts[0] ||
	insGrowCosts[0] != delGrowCosts[0]) {
      ERR("can't combine option -F > 0 with option -A or -B");
    }
  }
}

static int percent(int val, int percentage) {
  return (percentage == 100) ? val : val * percentage / 100;
}

void LastalArguments::setDefaultsFromMatrix( double lambda, int minScore ){
  if( outputType < 2 && minScoreGapped < 0 ) minScoreGapped = minScoreGapless;
  if( minScoreGapped < 0 ){
    if( outputType > 0 && minScore < 0 )
      ERR("\
can't calculate E-values: maybe the mismatch or gap costs are too weak.\n\
To proceed without E-values, set a score threshold with option -e.");
    minScoreGapped = minScore;
  }
  if( outputType < 2 ) minScoreGapless = minScoreGapped;

  if( temperature < 0 ) temperature = 1 / lambda;

  int defaultMaxScoreDrop = std::max(minScoreGapped - 1, 0);
  defaultMaxScoreDrop = std::max(defaultMaxScoreDrop, maxDropGapless);

  if (maxDropFinalSuffix == '%') {
    maxDropFinal = percent(defaultMaxScoreDrop, maxDropFinal);
  } else if (maxDropFinalSuffix == 'g') {
    maxDropFinal = std::min(minGapCost(maxDropFinal), defaultMaxScoreDrop);
  }

  if (maxDropGappedSuffix == '%') {
    maxDropGapped = percent(maxDropFinal, maxDropGapped);
  } else if (maxDropGappedSuffix == 'g') {
    maxDropGapped = std::min(minGapCost(maxDropGapped), maxDropFinal);
  }

  if( maxDropGapless < 0 ){  // should it depend on temperature or ...?
    if( temperature < 0 ) maxDropGapless = 0;  // shouldn't happen
    else                  maxDropGapless = int( 10.0 * temperature + 0.5 );
    maxDropGapless = std::min( maxDropGapless, maxDropGapped );
  }
}

void LastalArguments::writeCommented( std::ostream& stream ) const{
  stream << '#';
  stream << " a="; writeIntList(stream, delOpenCosts);
  stream << " b="; writeIntList(stream, delGrowCosts);
  stream << " A="; writeIntList(stream, insOpenCosts);
  stream << " B="; writeIntList(stream, insGrowCosts);
  if( gapPairCost > 0 )
    stream << " c=" << gapPairCost;
  if( isTranslated() )
    stream << " F=" << frameshiftCost;
  stream << " e=" << minScoreGapped;
  stream << " d=" << minScoreGapless;
  stream << " x=" << maxDropGapped;
  stream << " y=" << maxDropGapless;
  stream << " z=" << maxDropFinal;
  stream << " D=" << queryLettersPerRandomAlignment;
  stream << " E=" << maxEvalue;
  stream << '\n';

  stream << '#';
  stream << " R=" << isKeepLowercase << tantanSetting;
  stream << " u=" << maskLowercase;
  stream << " s=" << strand;
  stream << " S=" << isQueryStrandMatrix;
  stream << " M=" << isGreedy;
  stream << " T=" << globality;
  stream << " m=" << oneHitMultiplicity;
  stream << " l=" << minHitDepth;
  if( maxHitDepth + 1 > 0 )
    stream << " L=" << maxHitDepth;
  stream << " n=" << maxGaplessAlignmentsPerQueryPosition;
  if( maxAlignmentsPerQueryStrand + 1 > 0 )
    stream << " N=" << maxAlignmentsPerQueryStrand;
  if( cullingLimitForGaplessAlignments )
    stream << " C=" << cullingLimitForGaplessAlignments;
  if( cullingLimitForFinalAlignments )
    stream << " K=" << cullingLimitForFinalAlignments;
  stream << " k=" << queryStep;
  if( minimizerWindow > 1 )
    stream << " W=" << minimizerWindow;
  stream << " w=" << maxRepeatDistance;
  stream << " t=" << temperature;
  if( outputType > 4 && outputType < 7 )
    stream << " g=" << gamma;
  stream << " j=" << outputType;
  stream << " Q=" << (inputFormat % sequenceFormat::fasta);
  stream << '\n';

  stream << "# " << lastdbName << '\n';
}

}  // end namespace cbrc
