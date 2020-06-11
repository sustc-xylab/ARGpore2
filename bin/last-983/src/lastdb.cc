// Copyright 2008, 2009, 2010, 2011, 2013, 2014 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "last.hh"

#include "LastdbArguments.hh"
#include "SubsetSuffixArray.hh"
#include "TantanMasker.hh"
#include "zio.hh"
#include "stringify.hh"
#include "threadUtil.hh"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <numeric>  // accumulate

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << args.programName << ": " << x << '\n'

using namespace cbrc;

typedef unsigned long long countT;

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, const LastdbArguments& args ){
  if( !args.userAlphabet.empty() )  alph.fromString( args.userAlphabet );
  else if( args.isProtein )         alph.fromString( alph.protein );
  else                              alph.fromString( alph.dna );
}

// Does the first sequence look like it isn't really DNA?
bool isDubiousDna( const Alphabet& alph, const MultiSequence& multi ){
  const uchar* seq = multi.seqReader() + multi.seqBeg(0);
  unsigned dnaCount = 0;

  for( unsigned i = 0; i < 100; ++i ){  // look at the first 100 letters
    uchar c = alph.numbersToUppercase[ seq[i] ];
    if( c == alph.size ) return false;  // we hit the end of the sequence early
    if( c < alph.size || c == alph.encode[ (uchar)'N' ] ) ++dnaCount;
  }

  if( dnaCount < 90 ) return true;  // more than 10% unexpected letters
  else return false;
}

static void addSeeds( std::vector< CyclicSubsetSeed >& seeds,
		      const std::string& seedText,
		      const LastdbArguments& args, const Alphabet& alph ){
  std::istringstream iss( seedText );
  std::vector< std::string > seedAlphabet;
  std::string pattern;
  while( CyclicSubsetSeed::nextPattern( iss, seedAlphabet, pattern ) ){
    CyclicSubsetSeed s;
    s.init( seedAlphabet, pattern, args.isCaseSensitive, alph.encode );
    seeds.push_back(s);
  }
}

// Set up the seed pattern(s)
static void makeSubsetSeeds( std::vector< CyclicSubsetSeed >& seeds,
			     const std::string& seedText,
			     const LastdbArguments& args,
			     const Alphabet& alph ){
  const std::string& a = alph.letters;

  if( !args.subsetSeedFile.empty() ){
    addSeeds( seeds, seedText, args, alph );
  }
  else if( !args.seedPatterns.empty() ){
    for( unsigned x = 0; x < args.seedPatterns.size(); ++x ){
      const std::string& p = args.seedPatterns[x];
      std::string s = CyclicSubsetSeed::stringFromPatterns( p, a );
      addSeeds( seeds, s, args, alph );
    }
  }
  else{
    std::string s = (alph.letters == alph.dna)
      ? CyclicSubsetSeed::stringFromName( "YASS" )
      : CyclicSubsetSeed::stringFromPatterns( "1", a );
    addSeeds( seeds, s, args, alph );
  }

  if( seeds.empty() ) ERR( "no seed patterns" );
}

void writeLastalOptions( std::ostream& out, const std::string& seedText ){
  std::string trigger = "#lastal";
  std::istringstream iss( seedText );
  std::string line;
  while( getline( iss, line ) )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      out << line << '\n';
}

void writePrjFile( const std::string& fileName, const LastdbArguments& args,
		   const Alphabet& alph, countT sequenceCount,
		   const std::vector<countT>& letterCounts,
		   bool isFastq, unsigned volumes, unsigned numOfIndexes,
		   const std::string& seedText ){
  countT letterTotal = std::accumulate( letterCounts.begin(),
                                        letterCounts.end(), countT(0) );

  std::ofstream f( fileName.c_str() );
  f << "version=" <<
#include "version.hh"
    << '\n';
  f << "alphabet=" << alph << '\n';
  f << "numofsequences=" << sequenceCount << '\n';
  f << "numofletters=" << letterTotal << '\n';
  f << "letterfreqs=";
  for( unsigned i = 0; i < letterCounts.size(); ++i ){
    if( i > 0 ) f << ' ';
    f << letterCounts[i];
  }
  f << '\n';

  if( !args.isCountsOnly ){
    f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
    f << "keeplowercase=" << args.isKeepLowercase << '\n';
    if( args.tantanSetting ){
      f << "tantansetting=" << args.tantanSetting << '\n';
    }
    f << "masklowercase=" << args.isCaseSensitive << '\n';
    if( isFastq ){
      f << "sequenceformat=" << args.inputFormat << '\n';
    }
    if( args.minimizerWindow > 1 ){
      // Maybe this should be written (and read) by the indexes, so
      // each index can have a different window?
      f << "minimizerwindow=" << args.minimizerWindow << '\n';
    }
    if( volumes+1 > 0 ){
      f << "volumes=" << volumes << '\n';
    }
    else{
      f << "numofindexes=" << numOfIndexes << '\n';
    }
    f << "integersize=" << (sizeof(indexT) * CHAR_BIT) << '\n';
    writeLastalOptions( f, seedText );
  }

  f.close();
  if( !f ) ERR( "can't write file: " + fileName );
}

static void preprocessSomeSeqs(MultiSequence *multi,
			       const TantanMasker *masker,
			       const uchar *maskTable,
			       size_t numOfChunks,
			       size_t chunkNum) {
  size_t beg = firstSequenceInChunk(*multi, numOfChunks, chunkNum);
  size_t end = firstSequenceInChunk(*multi, numOfChunks, chunkNum + 1);
  uchar *w = multi->seqWriter();
  for (size_t i = beg; i < end; ++i)
    masker->mask(w + multi->seqBeg(i), w + multi->seqEnd(i), maskTable);
}

static void preprocessSeqs(MultiSequence &multi,
			   const TantanMasker &masker,
			   const uchar *maskTable,
			   size_t numOfChunks) {
#ifdef HAS_CXX_THREADS
  std::vector<std::thread> threads(numOfChunks - 1);
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1] = std::thread(preprocessSomeSeqs,
				 &multi, &masker, maskTable, numOfChunks, i);
#endif
  preprocessSomeSeqs(&multi, &masker, maskTable, numOfChunks, 0);
#ifdef HAS_CXX_THREADS
  for (size_t i = 1; i < numOfChunks; ++i)
    threads[i - 1].join();
#endif
}

// Make one database volume, from one batch of sequences
void makeVolume( std::vector< CyclicSubsetSeed >& seeds,
		 MultiSequence& multi, const LastdbArguments& args,
		 const Alphabet& alph, std::vector<countT>& letterCounts,
		 const TantanMasker& masker, unsigned numOfThreads,
		 const std::string& seedText, const std::string& baseName ){
  size_t numOfIndexes = seeds.size();
  size_t numOfSequences = multi.finishedSequences();
  size_t textLength = multi.seqBeg(numOfSequences);
  const uchar* seq = multi.seqReader();

  std::vector<countT> letterCountsInThisVolume(alph.size);
  alph.count(seq, seq + textLength, &letterCountsInThisVolume[0]);
  for (unsigned c = 0; c < alph.size; ++c) {
    letterCounts[c] += letterCountsInThisVolume[c];
  }

  if (args.isCountsOnly) return;

  if( args.tantanSetting ){
    LOG( "masking..." );
    preprocessSeqs( multi, masker, alph.numbersToLowercase, numOfThreads );
  }

  LOG( "writing..." );
  writePrjFile( baseName + ".prj", args, alph, numOfSequences,
		letterCountsInThisVolume, multi.qualsPerLetter(), -1,
		numOfIndexes, seedText );
  multi.toFiles( baseName );

  for( unsigned x = 0; x < numOfIndexes; ++x ){
    SubsetSuffixArray myIndex;
    seeds[x].swap( myIndex.getSeed() );

    LOG( "gathering..." );
    for( size_t i = 0; i < numOfSequences; ++i ){
      myIndex.addPositions( seq, multi.seqBeg(i), multi.seqEnd(i),
			    args.indexStep, args.minimizerWindow );
    }

    LOG( "sorting..." );
    myIndex.sortIndex( seq, args.minSeedLimit, args.childTableType,
		       numOfThreads );

    LOG( "bucketing..." );
    myIndex.makeBuckets( seq, args.bucketDepth );

    LOG( "writing..." );
    if( numOfIndexes > 1 ){
      myIndex.toFiles( baseName + char('a' + x), false, textLength );
    }
    else{
      myIndex.toFiles( baseName, true, textLength );
    }

    seeds[x].swap( myIndex.getSeed() );
  }

  LOG( "done!" );
}

// The max number of sequence letters, such that the total volume size
// is likely to be less than volumeSize bytes.  (This is crude, it
// neglects memory for the sequence names, and the fact that
// lowercase-masked letters and DNA "N"s aren't indexed.)
static indexT maxLettersPerVolume( const LastdbArguments& args,
				   size_t qualityCodesPerLetter,
				   unsigned numOfIndexes ){
  size_t bytesPerLetter = 1 + qualityCodesPerLetter;
  size_t maxIndexBytesPerPosition = sizeof(indexT) + 1;
  maxIndexBytesPerPosition *= numOfIndexes;
  size_t x = bytesPerLetter * args.indexStep + maxIndexBytesPerPosition;
  size_t y = args.volumeSize / x * args.indexStep;
  indexT z = y;
  if( z < y ) z = indexT(-1);
  return z;
}

static bool isRoomToDuplicateTheLastSequence(const MultiSequence &multi,
					     size_t maxSeqLen) {
  size_t n = multi.finishedSequences();
  size_t s = multi.seqBeg(n);
  return s <= maxSeqLen && s - multi.seqBeg(n - 1) <= maxSeqLen - s;
}

void lastdb( int argc, char** argv ){
  LastdbArguments args;
  args.fromArgs( argc, argv );

  std::string seedText;
  if( !args.subsetSeedFile.empty() ){
    seedText = CyclicSubsetSeed::stringFromName( args.subsetSeedFile );
    args.resetCumulativeOptions();
    args.fromString( seedText );  // read options from the seed file
    args.fromArgs( argc, argv );  // command line overrides seed file
  }

  unsigned numOfThreads =
    decideNumberOfThreads(args.numOfThreads, args.programName, args.verbosity);
  Alphabet alph;
  makeAlphabet( alph, args );
  TantanMasker tantanMasker;
  if( args.tantanSetting )
    tantanMasker.init( alph.isProtein(), args.tantanSetting > 1,
		       alph.letters, alph.encode );
  std::vector< CyclicSubsetSeed > seeds;
  makeSubsetSeeds( seeds, seedText, args, alph );
  MultiSequence multi;
  multi.initForAppending(1);
  alph.tr(multi.seqWriter(), multi.seqWriter() + multi.seqBeg(0));
  unsigned volumeNumber = 0;
  countT sequenceCount = 0;
  std::vector<countT> letterCounts( alph.size );
  indexT maxSeqLen = 0;

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;

  for( char** i = *inputBegin ? inputBegin : defaultInput; *i; ++i ){
    mcf::izstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while (appendSequence(multi, in, maxSeqLen, args.inputFormat, alph,
			  args.isKeepLowercase, 0)) {
      if (sequenceCount == 0) {
	maxSeqLen =
	  maxLettersPerVolume(args, multi.qualsPerLetter(), seeds.size());
	if (!args.isProtein && args.userAlphabet.empty() &&
	    isDubiousDna(alph, multi)) {
	  std::cerr << args.programName << ": that's some funny-lookin DNA\n";
	}
      }

      if( multi.isFinished() ){
	if (args.strand != 1) {
	  if (args.strand == 2) {
	    ++sequenceCount;
	    if (isRoomToDuplicateTheLastSequence(multi, maxSeqLen)) {
	      size_t lastSeq = multi.finishedSequences() - 1;
	      multi.duplicateOneSequence(lastSeq);
	    } else {
	      std::string baseName =
		args.lastdbName + stringify(volumeNumber++);
	      makeVolume(seeds, multi, args, alph, letterCounts,
			 tantanMasker, numOfThreads, seedText, baseName);
	      multi.eraseAllButTheLastSequence();
	    }
	  }
	  size_t lastSeq = multi.finishedSequences() - 1;
	  multi.reverseComplementOneSequence(lastSeq, alph.complement);
	}
        ++sequenceCount;
      }
      else{
	std::string baseName = args.lastdbName + stringify(volumeNumber++);
	makeVolume( seeds, multi, args, alph, letterCounts,
		    tantanMasker, numOfThreads, seedText, baseName );
	multi.reinitForAppending();
      }
    }
  }

  if( multi.finishedSequences() > 0 ){
    if( volumeNumber == 0 && !args.isCountsOnly ){
      makeVolume( seeds, multi, args, alph, letterCounts,
		  tantanMasker, numOfThreads, seedText, args.lastdbName );
      return;
    }
    std::string baseName = args.lastdbName + stringify(volumeNumber++);
    makeVolume( seeds, multi, args, alph, letterCounts,
		tantanMasker, numOfThreads, seedText, baseName );
  }

  writePrjFile( args.lastdbName + ".prj", args, alph, sequenceCount,
		letterCounts, multi.qualsPerLetter(), volumeNumber,
		seeds.size(), seedText );
}

int main( int argc, char** argv )
try{
  lastdb( argc, argv );
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
