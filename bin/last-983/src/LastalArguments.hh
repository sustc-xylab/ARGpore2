// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// This struct holds the command line arguments for lastal.

#ifndef LASTALARGUMENTS_HH
#define LASTALARGUMENTS_HH

#include "SequenceFormat.hh"

#include <climits>
#include <string>
#include <iosfwd>
#include <vector>
#include <stddef.h>  // size_t

namespace cbrc{

struct LastalArguments{
  // set the parameters to their default values:
  LastalArguments();

  // set parameters from a list of arguments:
  void fromArgs( int argc, char** argv, bool optionsOnly = false );

  // set parameters from a command line (by splitting it into arguments):
  void fromLine( const std::string& line );

  // set parameters from lines beginning with "#last":
  void fromString( const std::string& s );

  void resetCumulativeOptions() { verbosity = 0; }

  // get the name of the substitution score matrix:
  const char* matrixName( bool isProtein ) const;

  // set default option values that depend on input files:
  void setDefaultsFromAlphabet( bool isDna, bool isProtein,
				double numLettersInReference,
				bool isKeepRefLowercase, int refTantanSetting,
                                bool isCaseSensitiveSeeds, bool isVolumes,
				size_t refMinimizerWindow,
				unsigned realNumOfThreads );
  void setDefaultsFromMatrix( double lambda, int minScore );

  // write the parameter settings, starting each line with "#":
  void writeCommented( std::ostream& stream ) const;

  // are we doing translated alignment (DNA versus protein)?
  bool isTranslated() const{ return frameshiftCost >= 0; }

  // are we doing translated alignment with frameshifts?
  bool isFrameshift() const{ return frameshiftCost > 0; }

  // how many strands are we scanning (1 or 2)?
  int numOfStrands() const{ return (strand == 2) ? 2 : 1; }

  int minGapCost(int gapLength) const {
    int m = INT_MAX;
    for (size_t i = 0; i < delOpenCosts.size(); ++i) {
      m = std::min(m, delOpenCosts[i] + gapLength * delGrowCosts[i]);
    }
    for (size_t i = 0; i < insOpenCosts.size(); ++i) {
      m = std::min(m, insOpenCosts[i] + gapLength * insGrowCosts[i]);
    }
    return m;
  }

  // options:
  int outputFormat;
  int outputType;
  int strand;
  bool isQueryStrandMatrix;
  bool isGreedy;
  int globality;  // type of alignment: local, semi-global, etc.
  bool isKeepLowercase;
  int tantanSetting;
  int maskLowercase;
  double maxEvalue;
  double queryLettersPerRandomAlignment;
  int minScoreGapped;
  int minScoreGapless;
  int matchScore;
  int mismatchCost;
  std::vector<int> delOpenCosts;
  std::vector<int> delGrowCosts;
  std::vector<int> insOpenCosts;
  std::vector<int> insGrowCosts;
  int gapPairCost;
  int frameshiftCost;
  std::string matrixFile;
  int ambiguousLetterOpt;
  int maxDropGapped;
  char maxDropGappedSuffix;
  int maxDropGapless;
  int maxDropFinal;
  char maxDropFinalSuffix;
  sequenceFormat::Enum inputFormat;
  size_t minHitDepth;
  size_t maxHitDepth;
  size_t oneHitMultiplicity;
  size_t maxGaplessAlignmentsPerQueryPosition;
  size_t maxAlignmentsPerQueryStrand;
  size_t cullingLimitForGaplessAlignments;
  size_t cullingLimitForFinalAlignments;
  size_t queryStep;
  size_t minimizerWindow;
  size_t batchSize;  // approx size of query sequences to scan in 1 batch
  unsigned numOfThreads;
  size_t maxRepeatDistance;  // suppress repeats <= this distance apart
  double temperature;  // probability = exp( score / temperature ) / Z
  double gamma;        // parameter for gamma-centroid alignment
  std::string geneticCodeFile;
  int verbosity;

  // positional arguments:
  const char* programName;
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace cbrc
#endif  // LASTALARGUMENTS_HH
