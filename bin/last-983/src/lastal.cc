// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "last.hh"

#include "LastalArguments.hh"
#include "QualityPssmMaker.hh"
#include "OneQualityScoreMatrix.hh"
#include "TwoQualityScoreMatrix.hh"
#include "LastEvaluer.hh"
#include "GeneticCode.hh"
#include "SubsetMinimizerFinder.hh"
#include "SubsetSuffixArray.hh"
#include "Centroid.hh"
#include "AlignmentPot.hh"
#include "Alignment.hh"
#include "SegmentPairPot.hh"
#include "SegmentPair.hh"
#include "ScoreMatrix.hh"
#include "TantanMasker.hh"
#include "DiagonalTable.hh"
#include "GreedyXdropAligner.hh"
#include "gaplessXdrop.hh"
#include "gaplessPssmXdrop.hh"
#include "gaplessTwoQualityXdrop.hh"
#include "mcf_substitution_matrix_stats.hh"
#include "zio.hh"
#include "stringify.hh"
#include "threadUtil.hh"
#include <iomanip>  // setw
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << args.programName << ": " << x << '\n'
#define LOG2(x) if( args.verbosity > 1 ) std::cerr << args.programName << ": " << x << '\n'

static void warn( const char* programName, const char* s ){
  std::cerr << programName << ": " << s << '\n';
}

using namespace cbrc;

struct LastAligner {  // data that changes between queries
  Centroid centroid;
  GreedyXdropAligner greedyAligner;
  std::vector<int> qualityPssm;
  std::vector<AlignmentText> textAlns;
};

struct SubstitutionMatrices {
  int scores[scoreMatrixRowSize][scoreMatrixRowSize];
  int scoresMasked[scoreMatrixRowSize][scoreMatrixRowSize];
  mcf::SubstitutionMatrixStats stats;
  OneQualityScoreMatrix oneQual;
  OneQualityScoreMatrix oneQualMasked;
  OneQualityExpMatrix oneQualExp;
  QualityPssmMaker maker;
  TwoQualityScoreMatrix twoQual;
  TwoQualityScoreMatrix twoQualMasked;
};

namespace {
  typedef unsigned long long countT;

  LastalArguments args;
  Alphabet alph;
  Alphabet queryAlph;  // for translated alignment
  TantanMasker tantanMasker;
  GeneticCode geneticCode;
  const unsigned maxNumOfIndexes = 16;
  SubsetSuffixArray suffixArrays[maxNumOfIndexes];
  ScoreMatrix scoreMatrix;
  SubstitutionMatrices fwdMatrices;
  SubstitutionMatrices revMatrices;
  mcf::GapCosts gapCosts;
  std::vector<LastAligner> aligners;
  LastEvaluer evaluer;
  MultiSequence query;  // sequence that hasn't been indexed by lastdb
  MultiSequence text;  // sequence that has been indexed by lastdb
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  sequenceFormat::Enum referenceFormat = sequenceFormat::fasta;
  int minScoreGapless;
  int isCaseSensitiveSeeds = -1;  // initialize it to an "error" value
  unsigned numOfIndexes = 1;  // assume this value, if unspecified
}

void complementMatrix(const ScoreMatrixRow *from, ScoreMatrixRow *to) {
  for (unsigned i = 0; i < scoreMatrixRowSize; ++i)
    for (unsigned j = 0; j < scoreMatrixRowSize; ++j)
      to[i][j] = from[alph.complement[i]][alph.complement[j]];
}

// Meaningless for PSSMs, unless they have the same scale as the score matrix
static void
calculateSubstitutionScoreMatrixStatistics(const std::string &matrixName) {
  int *scoreMat[scoreMatrixRowSize];
  // the case-sensitivity of the matrix makes no difference here
  std::copy(scoreMatrix.caseSensitive,
	    scoreMatrix.caseSensitive + alph.size, scoreMat);

  mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
  if (args.temperature < 0) {
    const char *canonicalMatrixName = ScoreMatrix::canonicalName(matrixName);
    LOG("calculating matrix probabilities...");
    stats.calcUnbiased(canonicalMatrixName, scoreMat, alph.size);
    if (stats.isBad()) {
      static const char msg[] =
	"can't calculate probabilities: maybe the mismatch costs are too weak";
      if (isUseQuality(args.inputFormat) || args.outputType > 3) ERR(msg);
      LOG(msg);
      return;
    }
  } else {
    stats.calcFromScale(scoreMat, alph.size, args.temperature);
  }

  revMatrices.stats = stats;
  revMatrices.stats.flipDnaStrands(alph.complement);

  const double *p1 = stats.letterProbs1();
  const double *p2 = stats.letterProbs2();

  LOG( "matrix lambda=" << stats.lambda() );
  LOG( "matrix letter frequencies (upper=reference, lower=query):" );
  if( args.verbosity > 0 ){
    std::cerr << std::left;
    std::streamsize p = std::cerr.precision(2);
    unsigned e = alph.size;
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << alph.letters[i] << (i + 1 < e ? " " : "\n");
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << 100 * p1[i] << (i + 1 < e ? " " : "\n");
    for( unsigned i = 0; i < e; ++i )
      std::cerr << std::setw(3) << 100 * p2[i] << (i + 1 < e ? " " : "\n");
    std::cerr.precision(p);
    std::cerr << std::right;
  }
}

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( const std::string& matrixName,
		      const std::string& matrixFile ){
  if( !matrixName.empty() && !args.isGreedy ){
    scoreMatrix.fromString( matrixFile );
  }
  else{
    scoreMatrix.setMatchMismatch( args.matchScore, args.mismatchCost,
				  alph.letters );
  }

  scoreMatrix.init(alph.encode);
  if (args.outputType > 0) {
    calculateSubstitutionScoreMatrixStatistics(matrixName);
    const mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
    if (!stats.isBad()) {
      scoreMatrix.addAmbiguousScores(alph.letters == alph.dna,
				     args.ambiguousLetterOpt % 2,
				     args.ambiguousLetterOpt / 2,
				     alph.encode, stats.lambda(),
				     stats.letterProbs1(),
				     stats.letterProbs2());
      scoreMatrix.init(alph.encode);
    }
  }

  // If the input is a PSSM, the score matrix is not used, and its
  // maximum score should not be used.  Here, we try to set it to a
  // high enough value that it has no effect.  This is a kludge - it
  // would be nice to use the maximum PSSM score.
  if( args.inputFormat == sequenceFormat::pssm ) scoreMatrix.maxScore = 10000;
  // This would work, except the maxDrops aren't finalized yet:
  // maxScore = std::max(args.maxDropGapped, args.maxDropFinal) + 1;

  for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {  // xxx
    for (unsigned j = 0; j < scoreMatrixRowSize; ++j) {
      fwdMatrices.scores[i][j] = scoreMatrix.caseInsensitive[i][j];
      fwdMatrices.scoresMasked[i][j] = scoreMatrix.caseSensitive[i][j];
    }
  }

  if( args.isQueryStrandMatrix && args.strand != 1 ){
    complementMatrix(fwdMatrices.scores, revMatrices.scores);
    complementMatrix(fwdMatrices.scoresMasked, revMatrices.scoresMasked);
  }
}

void makeQualityScorers(SubstitutionMatrices &m, bool isCheck) {
  bool isMatchMismatch = (args.matrixFile.empty() && args.matchScore > 0);
  bool isPhred1 = isPhred( referenceFormat );
  int offset1 = qualityOffset( referenceFormat );
  bool isPhred2 = isPhred( args.inputFormat );
  int offset2 = qualityOffset( args.inputFormat );
  const uchar *toUnmasked = alph.numbersToUppercase;

  if( !isUseFastq( referenceFormat ) ){
    if( isUseFastq( args.inputFormat ) ){
      if ( isCheck ) LOG( "calculating per-quality scores..." );
      if( args.maskLowercase > 0 )
	m.oneQualMasked.init(m.scores, alph.size, m.stats,
			     isPhred2, offset2, toUnmasked, true);
      if( args.maskLowercase < 3 )
	m.oneQual.init(m.scores, alph.size, m.stats,
		       isPhred2, offset2, toUnmasked, false);
      const OneQualityScoreMatrix &q = (args.maskLowercase < 3) ?
	m.oneQual : m.oneQualMasked;
      if( args.outputType > 3 )
        m.oneQualExp.init(q, args.temperature);
      if( isCheck && args.verbosity > 0 )
	writeOneQualityScoreMatrix( q, alph.letters.c_str(),
				    offset2, std::cerr );
    }
    if( isUseQuality(args.inputFormat) ){
      m.maker.init(m.scores, alph.size, m.stats.lambda(),
		   isMatchMismatch, args.matchScore, -args.mismatchCost,
		   isPhred2, offset2, toUnmasked);
    }
  } else {
    if( isUseFastq( args.inputFormat ) ){
      if( args.maskLowercase > 0 )
	m.twoQualMasked.init(m.scores, m.stats,
			     isPhred1, offset1, isPhred2, offset2,
			     toUnmasked, true, isMatchMismatch);
      if( args.maskLowercase < 3 )
	m.twoQual.init(m.scores, m.stats,
		       isPhred1, offset1, isPhred2, offset2,
		       toUnmasked, false, isMatchMismatch);
      if( args.outputType > 3 )
        ERR( "fastq-versus-fastq column probabilities not implemented" );
    } else if (isCheck) {
      warn(args.programName,
	   "quality data not used for non-fastq query versus fastq reference");
    }
  }
}

void makeQualityScorers(){
  if( args.isGreedy ) return;

  if( args.isTranslated() )
    if( isUseQuality( args.inputFormat ) || isUseQuality( referenceFormat ) )
      return warn( args.programName,
		   "quality data not used for DNA-versus-protein alignment" );

  makeQualityScorers(fwdMatrices, true);
  if (args.isQueryStrandMatrix && args.strand != 1) {
    makeQualityScorers(revMatrices, false);
  }
}

// Calculate statistical parameters for the alignment scoring scheme
static void calculateScoreStatistics(const std::string& matrixName,
				     countT refLetters) {
  const mcf::SubstitutionMatrixStats &stats = fwdMatrices.stats;
  if (stats.isBad()) return;
  const char *canonicalMatrixName = ScoreMatrix::canonicalName( matrixName );
  if (args.temperature > 0 && !matrixName.empty()) canonicalMatrixName = " ";
  bool isGapped = (args.outputType > 1);
  bool isStandardGeneticCode = args.geneticCodeFile.empty();
  LOG( "getting E-value parameters..." );
  try{
    const mcf::GapCosts::Piece &del = gapCosts.delPieces[0];
    const mcf::GapCosts::Piece &ins = gapCosts.insPieces[0];
    evaluer.init( canonicalMatrixName, args.matchScore, args.mismatchCost,
                  alph.letters.c_str(), fwdMatrices.scoresMasked,
		  stats.letterProbs1(), stats.letterProbs2(), isGapped,
		  del.openCost, del.growCost, ins.openCost, ins.growCost,
                  args.frameshiftCost, geneticCode, isStandardGeneticCode,
		  args.verbosity );
    evaluer.setSearchSpace( refLetters, args.numOfStrands() );
    if( args.verbosity > 0 ) evaluer.writeParameters( std::cerr );
  }catch( const Sls::error& e ){
    LOG( "can't get E-value parameters for this scoring scheme" );
  }
}

// Read the .prj file for the whole database
void readOuterPrj( const std::string& fileName, unsigned& volumes,
                   size_t& refMinimizerWindow, size_t& minSeedLimit,
		   bool& isKeepRefLowercase, int& refTantanSetting,
                   countT& refSequences, countT& refLetters ){
  std::ifstream f( fileName.c_str() );
  if( !f ) ERR( "can't open file: " + fileName );
  unsigned version = 0;
  size_t fileBitsPerInt = 32;
  std::string trigger = "#lastal";

  std::string line, word;
  while( getline( f, line ) ){
    if( line.compare( 0, trigger.size(), trigger ) == 0 ){
      args.fromLine( line );
      continue;
    }
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "version" ) iss >> version;
    if( word == "alphabet" ) iss >> alph;
    if( word == "numofsequences" ) iss >> refSequences;
    if( word == "numofletters" ) iss >> refLetters;
    if( word == "maxunsortedinterval" ) iss >> minSeedLimit;
    if( word == "keeplowercase" ) iss >> isKeepRefLowercase;
    if( word == "tantansetting" ) iss >> refTantanSetting;
    if( word == "masklowercase" ) iss >> isCaseSensitiveSeeds;
    if( word == "sequenceformat" ) iss >> referenceFormat;
    if( word == "minimizerwindow" ) iss >> refMinimizerWindow;
    if( word == "volumes" ) iss >> volumes;
    if( word == "numofindexes" ) iss >> numOfIndexes;
    if( word == "integersize" ) iss >> fileBitsPerInt;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( alph.letters.empty() || refSequences+1 == 0 || refLetters+1 == 0 ||
      isCaseSensitiveSeeds < 0 || numOfIndexes > maxNumOfIndexes ||
      referenceFormat == sequenceFormat::prb ||
      referenceFormat == sequenceFormat::pssm ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) ERR( "can't read file: " + fileName );
  if( version < 294 && version > 0)
    ERR( "the lastdb files are old: please re-run lastdb" );

  if( fileBitsPerInt != sizeof(indexT) * CHAR_BIT ){
    if( fileBitsPerInt == 32 ) ERR( "please use lastal for " + fileName );
    if( fileBitsPerInt == 64 ) ERR( "please use lastal8 for " + fileName );
    ERR( "weird integersize in " + fileName );
  }
}

// Read a per-volume .prj file, with info about a database volume
void readInnerPrj( const std::string& fileName,
		   indexT& seqCount, indexT& seqLen ){
  std::ifstream f( fileName.c_str() );
  if( !f ) ERR( "can't open file: " + fileName );

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "numofsequences" ) iss >> seqCount;
    if( word == "numofletters" ) iss >> seqLen;
    if( word == "numofindexes" ) iss >> numOfIndexes;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( seqCount+1 == 0 || seqLen+1 == 0 || numOfIndexes > maxNumOfIndexes ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) ERR( "can't read file: " + fileName );
}

// Write match counts for each query sequence
void writeCounts( std::ostream& out ){
  for( indexT i = 0; i < matchCounts.size(); ++i ){
    out << query.seqName(i) << '\n';

    for( size_t j = args.minHitDepth; j < matchCounts[i].size(); ++j ){
      out << j << '\t' << matchCounts[i][j] << '\n';
    }

    out << '\n';  // blank line afterwards
  }
}

// Count all matches, of all sizes, of a query sequence against a suffix array
void countMatches( size_t queryNum, const uchar* querySeq ){
  size_t loopBeg = query.seqBeg(queryNum) - query.padBeg(queryNum);
  size_t loopEnd = query.seqEnd(queryNum) - query.padBeg(queryNum);
  if( args.minHitDepth > 1 )
    loopEnd -= std::min( args.minHitDepth - 1, loopEnd );

  for( size_t i = loopBeg; i < loopEnd; i += args.queryStep ){
    for( unsigned x = 0; x < numOfIndexes; ++x )
      suffixArrays[x].countMatches( matchCounts[queryNum], querySeq + i,
				    text.seqReader(), args.maxHitDepth );
  }
}

static const uchar *getQueryQual(size_t queryNum) {
  const uchar *q = query.qualityReader();
  if (q) q += query.padBeg(queryNum) * query.qualsPerLetter();
  return q;
}

static const ScoreMatrixRow *getQueryPssm(const LastAligner &aligner,
					  size_t queryNum) {
  if (args.isGreedy) return 0;
  if (args.inputFormat == sequenceFormat::pssm)
    return query.pssmReader() + query.padBeg(queryNum);
  const std::vector<int> &qualityPssm = aligner.qualityPssm;
  if (qualityPssm.empty())
    return 0;
  return reinterpret_cast<const ScoreMatrixRow *>(&qualityPssm[0]);
}

namespace Phase{ enum Enum{ gapless, gapped, final }; }

static bool isMaskLowercase(Phase::Enum e) {
  return (e < 1 && args.maskLowercase > 0) || args.maskLowercase > 2;
}

struct Dispatcher{
  const uchar* a;  // the reference sequence
  const uchar* b;  // the query sequence
  const uchar* i;  // the reference quality data
  const uchar* j;  // the query quality data
  const ScoreMatrixRow* p;  // the query PSSM
  const ScoreMatrixRow* m;  // the score matrix
  const TwoQualityScoreMatrix& t;
  int d;  // the maximum score drop
  int z;

  Dispatcher( Phase::Enum e, const LastAligner& aligner, size_t queryNum,
	      const SubstitutionMatrices &matrices, const uchar* querySeq ) :
      a( text.seqReader() ),
      b( querySeq ),
      i( text.qualityReader() ),
      j( getQueryQual(queryNum) ),
      p( getQueryPssm(aligner, queryNum) ),
      m( isMaskLowercase(e) ? matrices.scoresMasked : matrices.scores ),
      t( isMaskLowercase(e) ? matrices.twoQualMasked : matrices.twoQual ),
      d( (e == Phase::gapless) ? args.maxDropGapless :
         (e == Phase::gapped ) ? args.maxDropGapped : args.maxDropFinal ),
      z( t ? 2 : p ? 1 : 0 ){}

  int gaplessOverlap( indexT x, indexT y, size_t &rev, size_t &fwd ) const{
    if( z==0 ) return gaplessXdropOverlap( a+x, b+y, m, d, rev, fwd );
    if( z==1 ) return gaplessPssmXdropOverlap( a+x, p+y, d, rev, fwd );
    return gaplessTwoQualityXdropOverlap( a+x, i+x, b+y, j+y, t, d, rev, fwd );
  }

  int forwardGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return forwardGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return forwardGaplessPssmXdropScore( a+x, p+y, d );
    return forwardGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  int reverseGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return reverseGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return reverseGaplessPssmXdropScore( a+x, p+y, d );
    return reverseGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  indexT forwardGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return forwardGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return forwardGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return forwardGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  indexT reverseGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return reverseGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return reverseGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return reverseGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  bool isOptimalGapless( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return isOptimalGaplessXdrop( a+x, a+e, b+y, m, d );
    if( z==1 ) return isOptimalGaplessPssmXdrop( a+x, a+e, p+y, d );
    return isOptimalGaplessTwoQualityXdrop( a+x, a+e, i+x, b+y, j+y, t, d );
  }

  int gaplessScore( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return gaplessAlignmentScore( a+x, a+e, b+y, m );
    if( z==1 ) return gaplessPssmAlignmentScore( a+x, a+e, p+y );
    return gaplessTwoQualityAlignmentScore( a+x, a+e, i+x, b+y, j+y, t );
  }
};

static bool isCollatedAlignments() {
  return args.outputFormat == 'b' || args.outputFormat == 'B' ||
    args.cullingLimitForFinalAlignments;
}

static void printAndDelete(char *text) {
  std::cout << text;
  delete[] text;
}

static void writeAlignment(LastAligner &aligner, const Alignment &aln,
			   size_t queryNum, const uchar* querySeq,
			   const AlignmentExtras &extras = AlignmentExtras()) {
  AlignmentText a = aln.write(text, query, queryNum, querySeq,
			      args.isTranslated(), alph, evaluer,
			      args.outputFormat, extras);
  if (isCollatedAlignments() || aligners.size() > 1)
    aligner.textAlns.push_back(a);
  else
    printAndDelete(a.text);
}

static void writeSegmentPair(LastAligner &aligner, const SegmentPair &s,
			     size_t queryNum, const uchar* querySeq) {
  Alignment a;
  a.fromSegmentPair(s);
  writeAlignment(aligner, a, queryNum, querySeq);
}

// Find query matches to the suffix array, and do gapless extensions
void alignGapless( LastAligner& aligner, SegmentPairPot& gaplessAlns,
		   size_t queryNum, const SubstitutionMatrices &matrices,
		   const uchar* querySeq ){
  const bool isOverlap = (args.globality && args.outputType == 1);
  Dispatcher dis(Phase::gapless, aligner, queryNum, matrices, querySeq);
  DiagonalTable dt;  // record already-covered positions on each diagonal
  countT matchCount = 0, gaplessExtensionCount = 0, gaplessAlignmentCount = 0;
  size_t significantAlignmentCount = 0;

  size_t loopBeg = query.seqBeg(queryNum) - query.padBeg(queryNum);
  size_t loopEnd = query.seqEnd(queryNum) - query.padBeg(queryNum);
  if( args.minHitDepth > 1 )
    loopEnd -= std::min( args.minHitDepth - 1, loopEnd );

  std::vector< SubsetMinimizerFinder > minFinders( numOfIndexes );
  for( unsigned x = 0; x < numOfIndexes; ++x ){
    minFinders[x].init( suffixArrays[x].getSeed(), dis.b, loopBeg, loopEnd );
  }

  for( indexT i = loopBeg; i < loopEnd; i += args.queryStep ){
    for( unsigned x = 0; x < numOfIndexes; ++x ){
      const SubsetSuffixArray& sax = suffixArrays[x];
      if( args.minimizerWindow > 1 &&
	  !minFinders[x].isMinimizer( sax.getSeed(), dis.b, i, loopEnd,
				      args.minimizerWindow ) ) continue;
      const indexT* beg;
      const indexT* end;
      sax.match( beg, end, dis.b + i, dis.a,
		 args.oneHitMultiplicity, args.minHitDepth, args.maxHitDepth );
      matchCount += end - beg;

      // Tried: if we hit a delimiter when using contiguous seeds, then
      // increase "i" to the delimiter position.  This gave a speed-up
      // of only 3%, with 34-nt tags.

      size_t gaplessAlignmentsPerQueryPosition = 0;

      for( /* noop */; beg < end; ++beg ){  // loop over suffix-array matches
	if( gaplessAlignmentsPerQueryPosition ==
	    args.maxGaplessAlignmentsPerQueryPosition ) break;

	indexT j = *beg;  // coordinate in the reference sequence
	if( dt.isCovered( i, j ) ) continue;
	++gaplessExtensionCount;
	int score;

	if( isOverlap ){
	  size_t revLen, fwdLen;
	  score = dis.gaplessOverlap( j, i, revLen, fwdLen );
	  if( score < minScoreGapless ) continue;
	  SegmentPair sp( j - revLen, i - revLen, revLen + fwdLen, score );
	  dt.addEndpoint( sp.end2(), sp.end1() );
	  writeSegmentPair( aligner, sp, queryNum, querySeq );
	}else{
	  int fs = dis.forwardGaplessScore( j, i );
	  int rs = dis.reverseGaplessScore( j, i );
	  score = fs + rs;

	  // Tried checking the score after isOptimal & addEndpoint,
	  // but the number of extensions decreased by < 10%, and it
	  // was slower overall.
	  if( score < minScoreGapless ) continue;

	  indexT tEnd = dis.forwardGaplessEnd( j, i, fs );
	  indexT tBeg = dis.reverseGaplessEnd( j, i, rs );
	  indexT qBeg = i - (j - tBeg);
	  if( !dis.isOptimalGapless( tBeg, tEnd, qBeg ) ) continue;
	  SegmentPair sp( tBeg, qBeg, tEnd - tBeg, score );
	  dt.addEndpoint( sp.end2(), sp.end1() );

	  if( args.outputType == 1 ){  // we just want gapless alignments
	    writeSegmentPair( aligner, sp, queryNum, querySeq );
	  }
	  else{
	    gaplessAlns.add(sp);  // add the gapless alignment to the pot
	  }
	}

	++gaplessAlignmentsPerQueryPosition;
	++gaplessAlignmentCount;

	if( score >= args.minScoreGapped &&
	    ++significantAlignmentCount >= args.maxAlignmentsPerQueryStrand ) {
	  i = loopEnd;
	  break;
	}
      }
    }
  }

  LOG2( "initial matches=" << matchCount );
  LOG2( "gapless extensions=" << gaplessExtensionCount );
  LOG2( "gapless alignments=" << gaplessAlignmentCount );
}

// Shrink the SegmentPair to its longest run of identical matches.
// This trims off possibly unreliable parts of the gapless alignment.
// It may not be the best strategy for protein alignment with subset
// seeds: there could be few or no identical matches...
void shrinkToLongestIdenticalRun( SegmentPair& sp, const Dispatcher& dis ){
  sp.maxIdenticalRun( dis.a, dis.b, alph.numbersToUppercase );
  sp.score = dis.gaplessScore( sp.beg1(), sp.end1(), sp.beg2() );
}

// Do gapped extensions of the gapless alignments
void alignGapped( LastAligner& aligner,
		  AlignmentPot& gappedAlns, SegmentPairPot& gaplessAlns,
                  size_t queryNum, const SubstitutionMatrices &matrices,
		  const uchar* querySeq, Phase::Enum phase ){
  Dispatcher dis(phase, aligner, queryNum, matrices, querySeq);
  size_t frameSize = args.isFrameshift() ? (query.padLen(queryNum) / 3) : 0;
  countT gappedExtensionCount = 0, gappedAlignmentCount = 0;

  // Redo the gapless extensions, using gapped score parameters.
  // Without this, if we self-compare a huge sequence, we risk getting
  // huge gapped extensions.
  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.items[i];

    int fs = dis.forwardGaplessScore( sp.beg1(), sp.beg2() );
    int rs = dis.reverseGaplessScore( sp.beg1(), sp.beg2() );
    indexT tEnd = dis.forwardGaplessEnd( sp.beg1(), sp.beg2(), fs );
    indexT tBeg = dis.reverseGaplessEnd( sp.beg1(), sp.beg2(), rs );
    indexT qBeg = sp.beg2() - (sp.beg1() - tBeg);
    sp = SegmentPair( tBeg, qBeg, tEnd - tBeg, fs + rs );

    if( !dis.isOptimalGapless( tBeg, tEnd, qBeg ) ){
      SegmentPairPot::mark(sp);
    }
  }

  erase_if( gaplessAlns.items, SegmentPairPot::isMarked );

  gaplessAlns.cull( args.cullingLimitForGaplessAlignments );
  gaplessAlns.sort();  // sort by score descending, and remove duplicates

  LOG2( "redone gapless alignments=" << gaplessAlns.size() );

  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.get(i);

    if( SegmentPairPot::isMarked(sp) ) continue;

    Alignment aln;
    AlignmentExtras extras;  // not used
    aln.seed = sp;

    shrinkToLongestIdenticalRun( aln.seed, dis );

    // do gapped extension from each end of the seed:
    aln.makeXdrop( aligner.centroid, aligner.greedyAligner, args.isGreedy,
		   dis.a, dis.b, args.globality, dis.m, scoreMatrix.maxScore,
		   gapCosts, dis.d, args.frameshiftCost, frameSize,
		   dis.p, dis.t, dis.i, dis.j, alph, extras );
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if( !aln.isOptimal( dis.a, dis.b, args.globality, dis.m, dis.d, gapCosts,
			args.frameshiftCost, frameSize, dis.p,
                        dis.t, dis.i, dis.j ) ){
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-redundantization.
      continue;
    }

    gaplessAlns.markAllOverlaps( aln.blocks );
    gaplessAlns.markTandemRepeats( aln.seed, args.maxRepeatDistance );

    if( phase == Phase::final ) gappedAlns.add(aln);
    else SegmentPairPot::markAsGood(sp);

    ++gappedAlignmentCount;
    if( gappedAlignmentCount >= args.maxAlignmentsPerQueryStrand ) break;
  }

  LOG2( "gapped extensions=" << gappedExtensionCount );
  LOG2( "gapped alignments=" << gappedAlignmentCount );
}

// Print the gapped alignments, after optionally calculating match
// probabilities and re-aligning using the gamma-centroid algorithm
void alignFinish( LastAligner& aligner, const AlignmentPot& gappedAlns,
		  size_t queryNum, const SubstitutionMatrices &matrices,
		  const uchar* querySeq ){
  Centroid& centroid = aligner.centroid;
  Dispatcher dis(Phase::final, aligner, queryNum, matrices, querySeq);
  size_t queryLen = query.padLen(queryNum);
  size_t frameSize = args.isFrameshift() ? (queryLen / 3) : 0;

  if( args.outputType > 3 ){
    if( dis.p ){
      if (args.outputType == 7) {
	centroid.setScoreMatrix(dis.m, args.temperature);
	centroid.setLetterProbsPerPosition(alph.size, queryLen, dis.b, dis.j,
					   isUseFastq(args.inputFormat),
					   matrices.maker.qualToProbRight(),
					   matrices.stats.letterProbs2(),
					   alph.numbersToUppercase);
      }
      centroid.setPssm(dis.p, queryLen, args.temperature,
		       matrices.oneQualExp, dis.b, dis.j);
    }
    else{
      centroid.setScoreMatrix( dis.m, args.temperature );
    }
    centroid.setOutputType( args.outputType );
  }

  for( size_t i = 0; i < gappedAlns.size(); ++i ){
    const Alignment& aln = gappedAlns.items[i];
    if( args.outputType < 4 ){
      writeAlignment( aligner, aln, queryNum, querySeq );
    }
    else{  // calculate match probabilities:
      Alignment probAln;
      AlignmentExtras extras;
      probAln.seed = aln.seed;
      probAln.makeXdrop( centroid, aligner.greedyAligner, args.isGreedy,
			 dis.a, dis.b, args.globality,
			 dis.m, scoreMatrix.maxScore, gapCosts, dis.d,
                         args.frameshiftCost, frameSize,
			 dis.p, dis.t, dis.i, dis.j, alph, extras,
			 args.gamma, args.outputType );
      assert( aln.score != -INF );
      writeAlignment( aligner, probAln, queryNum, querySeq, extras );
    }
  }
}

static void eraseWeakAlignments(LastAligner &aligner, AlignmentPot &gappedAlns,
				size_t queryNum,
				const SubstitutionMatrices &matrices,
				const uchar *querySeq) {
  size_t frameSize = args.isFrameshift() ? (query.padLen(queryNum) / 3) : 0;
  Dispatcher dis(Phase::gapless, aligner, queryNum, matrices, querySeq);
  for (size_t i = 0; i < gappedAlns.size(); ++i) {
    Alignment &a = gappedAlns.items[i];
    if (!a.hasGoodSegment(dis.a, dis.b, args.minScoreGapped, dis.m, gapCosts,
			  args.frameshiftCost, frameSize,
			  dis.p, dis.t, dis.i, dis.j)) {
      AlignmentPot::mark(a);
    }
  }
  erase_if(gappedAlns.items, AlignmentPot::isMarked);
}

static bool lessForCulling(const AlignmentText &x, const AlignmentText &y) {
  if (x.strandNum != y.strandNum) return x.strandNum < y.strandNum;
  if (x.queryBeg  != y.queryBeg ) return x.queryBeg  < y.queryBeg;
  else                            return x.score     > y.score;
}

// Remove any alignment whose query range lies in LIMIT or more other
// alignments with higher score (and on the same strand):
static void cullFinalAlignments(std::vector<AlignmentText> &textAlns,
				size_t start) {
  if (!args.cullingLimitForFinalAlignments) return;
  sort(textAlns.begin() + start, textAlns.end(), lessForCulling);
  std::vector<size_t> stash;  // alignments that might dominate subsequent ones
  size_t i = start;  // number of kept alignments so far
  for (size_t j = start; j < textAlns.size(); ++j) {
    AlignmentText &x = textAlns[j];
    size_t numOfDominators = 0;  // number of alignments that dominate x
    size_t a = 0;  // number of kept stash-items so far
    for (size_t b = 0; b < stash.size(); ++b) {
      size_t k = stash[b];
      AlignmentText &y = textAlns[k];
      if (y.strandNum < x.strandNum) break;  // drop the stash
      if (y.queryEnd <= x.queryBeg) continue;  // drop this stash-item
      stash[a++] = k;  // keep this stash-item
      if (y.queryEnd >= x.queryEnd && y.score > x.score) ++numOfDominators;
    }
    stash.resize(a);
    if (numOfDominators >= args.cullingLimitForFinalAlignments) {
      delete[] x.text;
    } else {
      stash.push_back(i);
      textAlns[i++] = x;  // keep this alignment
    }
  }
  textAlns.resize(i);
}

static void printAndClear(std::vector<AlignmentText> &textAlns) {
  for (size_t i = 0; i < textAlns.size(); ++i)
    printAndDelete(textAlns[i].text);
  textAlns.clear();
}

static void printAndClearAll() {
  for (size_t i = 0; i < aligners.size(); ++i)
    printAndClear(aligners[i].textAlns);
}

void makeQualityPssm( LastAligner& aligner,
		      size_t queryNum, const SubstitutionMatrices &matrices,
		      const uchar* querySeq, bool isMask ){
  if (!isUseQuality(args.inputFormat) || isUseQuality(referenceFormat)) return;
  if( args.isTranslated() || args.isGreedy ) return;

  std::vector<int> &qualityPssm = aligner.qualityPssm;
  size_t queryLen = query.padLen(queryNum);
  qualityPssm.resize(queryLen * scoreMatrixRowSize);

  const uchar *seqBeg = querySeq;
  const uchar *seqEnd = seqBeg + queryLen;
  const uchar *q = getQueryQual(queryNum);
  int *pssm = &qualityPssm[0];

  if( args.inputFormat == sequenceFormat::prb ){
    matrices.maker.make( seqBeg, seqEnd, q, pssm, isMask );
  }
  else {
    const OneQualityScoreMatrix &m =
      isMask ? matrices.oneQualMasked : matrices.oneQual;
    makePositionSpecificScoreMatrix( m, seqBeg, seqEnd, q, pssm );
  }
}

// Scan one query sequence against one database volume
void scan(LastAligner& aligner, size_t queryNum,
	  const SubstitutionMatrices &matrices, const uchar* querySeq) {
  if( args.outputType == 0 ){  // we just want match counts
    countMatches( queryNum, querySeq );
    return;
  }

  bool isMask = (args.maskLowercase > 0);
  makeQualityPssm(aligner, queryNum, matrices, querySeq, isMask);

  SegmentPairPot gaplessAlns;
  alignGapless(aligner, gaplessAlns, queryNum, matrices, querySeq);
  if( args.outputType == 1 ) return;  // we just want gapless alignments
  if( gaplessAlns.size() == 0 ) return;

  if( args.maskLowercase == 1 || args.maskLowercase == 2 )
    makeQualityPssm(aligner, queryNum, matrices, querySeq, false);

  AlignmentPot gappedAlns;

  if( args.maxDropFinal != args.maxDropGapped ){
    alignGapped( aligner, gappedAlns, gaplessAlns,
		 queryNum, matrices, querySeq, Phase::gapped );
    erase_if( gaplessAlns.items, SegmentPairPot::isNotMarkedAsGood );
  }

  alignGapped( aligner, gappedAlns, gaplessAlns,
	       queryNum, matrices, querySeq, Phase::final );
  if( gappedAlns.size() == 0 ) return;

  if (args.maskLowercase == 2) {
    makeQualityPssm(aligner, queryNum, matrices, querySeq, true);
    eraseWeakAlignments(aligner, gappedAlns, queryNum, matrices, querySeq);
    LOG2("lowercase-filtered alignments=" << gappedAlns.size());
    if (gappedAlns.size() == 0) return;
    if (args.outputType > 3)
      makeQualityPssm(aligner, queryNum, matrices, querySeq, false);
  }

  if( args.outputType > 2 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG2( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  if( !isCollatedAlignments() ) gappedAlns.sort();  // sort by score
  alignFinish(aligner, gappedAlns, queryNum, matrices, querySeq);
}

static void tantanMaskOneQuery(size_t queryNum, uchar *querySeq) {
  size_t beg = query.seqBeg(queryNum) - query.padBeg(queryNum);
  size_t end = query.seqEnd(queryNum) - query.padBeg(queryNum);
  tantanMasker.mask(querySeq + beg, querySeq + end, alph.numbersToLowercase);
}

static void tantanMaskTranslatedQuery(size_t queryNum, uchar *querySeq) {
  size_t frameSize = query.padLen(queryNum) / 3;
  size_t dnaBeg = query.seqBeg(queryNum) - query.padBeg(queryNum);
  size_t dnaLen = query.seqLen(queryNum);
  for (int frame = 0; frame < 3; ++frame) {
    if (dnaLen < 3) break;
    size_t aaBeg = dnaToAa(dnaBeg++, frameSize);
    size_t aaLen = dnaLen-- / 3;
    size_t aaEnd = aaBeg + aaLen;
    tantanMasker.mask(querySeq + aaBeg, querySeq + aaEnd,
		      alph.numbersToLowercase);
  }
}

// Scan one query sequence strand against one database volume,
// after optionally translating and/or masking the query
void translateAndScan(LastAligner& aligner,
		      size_t queryNum, const SubstitutionMatrices &matrices) {
  const uchar* querySeq = query.seqReader() + query.padBeg(queryNum);
  std::vector<uchar> modifiedQuery;
  size_t size = query.padLen(queryNum);

  if( args.isTranslated() ){
    modifiedQuery.resize( size );
    geneticCode.translate( querySeq, querySeq + size, &modifiedQuery[0] );
    if( args.tantanSetting ){
      tantanMaskTranslatedQuery( queryNum, &modifiedQuery[0] );
    }
    querySeq = &modifiedQuery[0];
  }else{
    if( args.tantanSetting ){
      modifiedQuery.assign( querySeq, querySeq + size );
      tantanMaskOneQuery( queryNum, &modifiedQuery[0] );
      querySeq = &modifiedQuery[0];
    }
  }

  size_t oldNumOfAlns = aligner.textAlns.size();
  scan( aligner, queryNum, matrices, querySeq );
  cullFinalAlignments( aligner.textAlns, oldNumOfAlns );
}

static void alignOneQuery(LastAligner &aligner,
			  size_t queryNum, bool isFirstVolume) {
  if (args.strand == 2 && !isFirstVolume)
    query.reverseComplementOneSequence(queryNum, queryAlph.complement);

  if (args.strand != 0)
    translateAndScan(aligner, queryNum, fwdMatrices);

  if (args.strand == 2 || (args.strand == 0 && isFirstVolume))
    query.reverseComplementOneSequence(queryNum, queryAlph.complement);

  if (args.strand != 1)
    translateAndScan(aligner, queryNum,
		     args.isQueryStrandMatrix ? revMatrices : fwdMatrices);
}

static void alignSomeQueries(size_t chunkNum,
			     unsigned volume, unsigned volumeCount) {
  size_t numOfChunks = aligners.size();
  LastAligner &aligner = aligners[chunkNum];
  std::vector<AlignmentText> &textAlns = aligner.textAlns;
  size_t beg = firstSequenceInChunk(query, numOfChunks, chunkNum);
  size_t end = firstSequenceInChunk(query, numOfChunks, chunkNum + 1);
  bool isMultiVolume = (volumeCount > 1);
  bool isFirstVolume = (volume == 0);
  bool isFinalVolume = (volume + 1 == volumeCount);
  bool isFirstThread = (chunkNum == 0);
  bool isSort = isCollatedAlignments();
  bool isSortPerQuery = (isSort && !isMultiVolume);
  bool isPrintPerQuery = (isFirstThread && !(isSort && isMultiVolume));
  for (size_t i = beg; i < end; ++i) {
    size_t oldNumOfAlns = textAlns.size();
    alignOneQuery(aligner, i, isFirstVolume);
    if (isSortPerQuery) sort(textAlns.begin() + oldNumOfAlns, textAlns.end());
    if (isPrintPerQuery) printAndClear(textAlns);
  }
  if (isFinalVolume && isMultiVolume) {
    cullFinalAlignments(textAlns, 0);
    if (isSort) sort(textAlns.begin(), textAlns.end());
    if (isFirstThread) printAndClear(textAlns);
  }
}

static void scanOneVolume(unsigned volume, unsigned volumeCount) {
#ifdef HAS_CXX_THREADS
  size_t numOfChunks = aligners.size();
  std::vector<std::thread> threads(numOfChunks - 1);
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1] = std::thread(alignSomeQueries, i, volume, volumeCount);
  // Exceptions from threads are not handled nicely, but I don't
  // think it matters much.
#endif
  alignSomeQueries(0, volume, volumeCount);
#ifdef HAS_CXX_THREADS
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1].join();
#endif
}

void readIndex( const std::string& baseName, indexT seqCount ) {
  LOG( "reading " << baseName << "..." );
  text.fromFiles(baseName, seqCount, referenceFormat != sequenceFormat::fasta);
  for( unsigned x = 0; x < numOfIndexes; ++x ){
    if( numOfIndexes > 1 ){
      suffixArrays[x].fromFiles( baseName + char('a' + x),
				 isCaseSensitiveSeeds, alph.encode );
    }else{
      suffixArrays[x].fromFiles( baseName, isCaseSensitiveSeeds, alph.encode );
    }
  }
}

int calcMinScoreGapless(double numLettersInReference) {
  if (args.minScoreGapless >= 0) return args.minScoreGapless;

  // ***** Default setting for minScoreGapless *****

  // This attempts to ensure that the gapped alignment phase will be
  // reasonably fast relative to the gapless alignment phase.

  // The expected number of gapped extensions per query position is:
  // kGapless * referenceSize * exp(-lambdaGapless * minScoreGapless).

  // The number of gapless extensions per query position is
  // proportional to: maxGaplessAlignmentsPerQueryPosition * numOfIndexes.

  // So we want exp(lambdaGapless * minScoreGapless) to be
  // proportional to: kGapless * referenceSize / (n * numOfIndexes).

  // But we crudely ignore kGapless.

  // The proportionality constant was guesstimated by some limited
  // trial-and-error.  It should depend on the relative speeds of
  // gapless and gapped extensions.

  if (args.temperature < 0) return args.minScoreGapped;  // shouldn't happen

  double n = args.maxGaplessAlignmentsPerQueryPosition;
  if (args.maxGaplessAlignmentsPerQueryPosition + 1 == 0) n = 10;  // ?
  double x = 1000.0 * numLettersInReference / (n * numOfIndexes);
  if (x < 1) x = 1;
  int s = int(args.temperature * std::log(x) + 0.5);
  return std::min(s, args.minScoreGapped);
}

// Read one database volume
void readVolume( unsigned volumeNumber ){
  std::string baseName = args.lastdbName + stringify(volumeNumber);
  indexT seqCount = indexT(-1);
  indexT seqLen = indexT(-1);
  readInnerPrj( baseName + ".prj", seqCount, seqLen );
  minScoreGapless = calcMinScoreGapless(seqLen);
  readIndex( baseName, seqCount );
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes( unsigned volumes, std::ostream& out ){
  if( args.outputType == 0 ){
    matchCounts.clear();
    matchCounts.resize( query.finishedSequences() );
  }

  if( volumes+1 == 0 ) volumes = 1;
  bool isMultiVolume = (volumes > 1);

  for( unsigned i = 0; i < volumes; ++i ){
    if( text.unfinishedSize() == 0 || isMultiVolume ) readVolume( i );
    scanOneVolume( i, volumes );
    if( !isCollatedAlignments() ) printAndClearAll();
  }

  if( args.outputType == 0 ) writeCounts( out );
  printAndClearAll();
}

void writeHeader( countT refSequences, countT refLetters, std::ostream& out ){
  out << "# LAST version " <<
#include "version.hh"
      << "\n";
  out << "#\n";
  args.writeCommented( out );
  out << "# Reference sequences=" << refSequences
      << " normal letters=" << refLetters << "\n";
  if( args.outputType > 0 ) evaluer.writeCommented( out );
  out << "#\n";

  if( args.outputType == 0 ){  // we just want hit counts
    out << "# length\tcount\n"
	<< "#\n";
  }
  else{  // we want alignments
    if( args.inputFormat != sequenceFormat::pssm || !args.matrixFile.empty() ){
      // we're not reading PSSMs, or we bothered to specify a matrix file
      scoreMatrix.writeCommented( out );
      // Write lambda?
      out << "#\n";
    }

    if( args.outputFormat != 'b' && args.outputFormat != 'B' ) {
      out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
      out << "# in the reverse complement of the 2nd sequence are used.\n";
      out << "#\n";
    }

    if( args.outputFormat == 't' ){
      out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
	  << "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
    }
    if( args.outputFormat == 'm' ){
      out << "# name start alnSize strand seqSize alignment\n"
	  << "#\n";
    }
    if( args.outputFormat == 'b' || args.outputFormat == 'B' ){
      out << "# Fields: query id, subject id, % identity, alignment length, "
	  << "mismatches, gap opens, q. start, q. end, s. start, s. end";
      if( evaluer.isGood() ) out << ", evalue, bit score";
      if( args.outputFormat == 'B' ) out << ", query length, subject length";
      out << '\n';
    }
  }
}

void lastal( int argc, char** argv ){
  args.fromArgs( argc, argv );
  args.resetCumulativeOptions();  // because we will do fromArgs again

  unsigned volumes = unsigned(-1);
  size_t refMinimizerWindow = 1;  // assume this value, if not specified
  size_t minSeedLimit = 0;
  countT refSequences = -1;
  countT refLetters = -1;
  bool isKeepRefLowercase = true;
  int refTantanSetting = 0;
  readOuterPrj( args.lastdbName + ".prj", volumes,
		refMinimizerWindow, minSeedLimit,
		isKeepRefLowercase, refTantanSetting,
		refSequences, refLetters );
  bool isDna = (alph.letters == alph.dna);
  bool isProtein = alph.isProtein();

  args.fromArgs( argc, argv );  // command line overrides prj file

  std::string matrixName = args.matrixName( isProtein );
  std::string matrixFile;
  if( !matrixName.empty() ){
    matrixFile = ScoreMatrix::stringFromName( matrixName );
    args.resetCumulativeOptions();
    args.fromString( matrixFile );  // read options from the matrix file
    args.fromArgs( argc, argv );  // command line overrides matrix file
  }

  if( minSeedLimit > 1 ){
    if( args.outputType == 0 )
      ERR( "can't use option -j 0: need to re-run lastdb with i <= 1" );
    if( minSeedLimit > args.oneHitMultiplicity )
      ERR( "can't use option -m < " + stringify(minSeedLimit) +
	   ": need to re-run lastdb with i <= " +
	   stringify(args.oneHitMultiplicity) );
    if( args.minHitDepth > 1 )
      ERR( "can't use option -l > 1: need to re-run lastdb with i <= 1" );
  }

  aligners.resize( decideNumberOfThreads( args.numOfThreads,
					  args.programName, args.verbosity ) );
  bool isMultiVolume = (volumes+1 > 0 && volumes > 1);
  args.setDefaultsFromAlphabet( isDna, isProtein, refLetters,
				isKeepRefLowercase, refTantanSetting,
                                isCaseSensitiveSeeds, isMultiVolume,
				refMinimizerWindow, aligners.size() );
  if( args.tantanSetting )
    tantanMasker.init( isProtein, args.tantanSetting > 1,
		       alph.letters, alph.encode );
  makeScoreMatrix( matrixName, matrixFile );
  gapCosts.assign(args.delOpenCosts, args.delGrowCosts,
		  args.insOpenCosts, args.insGrowCosts, args.gapPairCost);

  if( args.isTranslated() ){
    if( isDna )  // allow user-defined alphabet
      ERR( "expected protein database, but got DNA" );
    queryAlph.fromString( queryAlph.dna );
    if( args.geneticCodeFile.empty() )
      geneticCode.fromString( geneticCode.standard );
    else
      geneticCode.fromFile( args.geneticCodeFile );
    geneticCode.codeTableSet( alph, queryAlph );
    query.initForAppending(3);
  }
  else{
    queryAlph = alph;
    query.initForAppending(1);
  }

  if( args.outputType > 0 ) calculateScoreStatistics( matrixName, refLetters );
  int minScore = evaluer.minScore( args.maxEvalue, 1e18 );
  args.setDefaultsFromMatrix( fwdMatrices.stats.lambda(), minScore );
  minScoreGapless = calcMinScoreGapless(refLetters);
  if( !isMultiVolume ) args.minScoreGapless = minScoreGapless;
  if( args.outputType > 0 ) makeQualityScorers();

  queryAlph.tr(query.seqWriter(), query.seqWriter() + query.seqBeg(0));

  if( volumes+1 == 0 ) readIndex( args.lastdbName, refSequences );

  std::ostream& out = std::cout;
  writeHeader( refSequences, refLetters, out );
  out.precision(3);  // print non-integers more compactly
  countT queryBatchCount = 0;
  countT sequenceCount = 0;
  indexT maxSeqLen = args.batchSize;
  if (maxSeqLen < args.batchSize) maxSeqLen = -1;

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;

  for( char** i = *inputBegin ? inputBegin : defaultInput; *i; ++i ){
    mcf::izstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while (appendSequence(query, in, maxSeqLen, args.inputFormat, queryAlph,
			  args.isKeepLowercase, args.maskLowercase > 1)) {
      if( query.isFinished() ){
	++sequenceCount;
      }else{
        // this enables downstream parsers to read one batch at a time:
        out << "# batch " << queryBatchCount++ << "\n";
	scanAllVolumes( volumes, out );
	query.reinitForAppending();
      }
    }
  }

  if( query.finishedSequences() > 0 ){
    out << "# batch " << queryBatchCount << "\n";
    scanAllVolumes( volumes, out );
  }

  out << "# Query sequences=" << sequenceCount << "\n";
}

int main( int argc, char** argv )
try{
  lastal( argc, argv );
  if (!flush(std::cout)) ERR( "write error" );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << argv[0] << ": out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << argv[0] << ": " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
