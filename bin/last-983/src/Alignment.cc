// Copyright 2008, 2009, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "Alphabet.hh"
#include "Centroid.hh"
#include "GeneticCode.hh"
#include "GreedyXdropAligner.hh"
#include "TwoQualityScoreMatrix.hh"
#include <cassert>

// make C++ tolerable:
#define IT(type) std::vector<type>::iterator

using namespace cbrc;

void Alignment::fromSegmentPair( const SegmentPair& sp ){
  blocks.assign( 1, sp );
  score = sp.score;
}

static void addExpectedCounts( double* expectedCounts,
			       const ExpectedCount& ec,
			       const Alphabet& alph ){
  for( unsigned i = 0; i < scoreMatrixRowSize; ++i ){
    unsigned x = alph.numbersToUppercase[i];
    if( x >= alph.size ) continue;
    for( unsigned j = 0; j < scoreMatrixRowSize; ++j ){
      unsigned y = alph.numbersToUppercase[j];
      if( y >= alph.size ) continue;
      expectedCounts[ x * alph.size + y ] += ec.emit[i][j];
    }
  }

  const int numEmissionCounts = alph.size * alph.size;
  double* transitionCounts = &expectedCounts[ numEmissionCounts ];

  transitionCounts[0] += ec.toMatch;
  transitionCounts[1] += ec.DD + ec.MD + ec.PD;  // deleted letter count
  transitionCounts[2] += ec.II + ec.MI + ec.DI + ec.PI;  // ins. letter count
  transitionCounts[3] += ec.MD + ec.PD;  // deletion open/close count
  transitionCounts[4] += ec.MI + ec.DI + ec.PI;  // insertion open/close count
  transitionCounts[5] += ec.DI;  // adjacent insertion & deletion count
  transitionCounts[7] += ec.PP + ec.MP;  // unaligned letter pair count
  transitionCounts[6] += ec.MP;  // pair-gap open/close count
  transitionCounts[8] += ec.PD;
  transitionCounts[9] += ec.PI;
  // MD = DM + DI - PD + DQ
  // MI = IM - DI - PI + IQ
  // PM = MP - PD - PI - PQ
  // DM + IM + PM = MD + MI + MP - DQ - IQ - PQ
}

static void countSeedMatches( double* expectedCounts,
			      const uchar* seq1beg, const uchar* seq1end,
			      const uchar* seq2beg, const Alphabet& alph ){
  while( seq1beg < seq1end ){
    unsigned x1 = alph.numbersToUppercase[ *seq1beg++ ];
    unsigned x2 = alph.numbersToUppercase[ *seq2beg++ ];
    if( x1 < alph.size && x2 < alph.size )
      ++expectedCounts[ x1 * alph.size + x2 ];
  }
}

// Does x precede and touch y in both sequences?
static bool isNext( const SegmentPair& x, const SegmentPair& y ){
  return x.end1() == y.beg1() && x.end2() == y.beg2();
}

void Alignment::makeXdrop( Centroid& centroid,
			   GreedyXdropAligner& greedyAligner, bool isGreedy,
			   const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix, int smMax,
			   const mcf::GapCosts& gap, int maxDrop,
			   int frameshiftCost, size_t frameSize,
			   const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2,
			   const Alphabet& alph, AlignmentExtras& extras,
			   double gamma, int outputType ){
  score = seed.score;
  if( outputType > 3 ) extras.fullScore = seed.score;

  if( outputType == 7 ){
    assert( seed.size > 0 );  // makes things easier to understand
    const int numEmissionCounts = alph.size * alph.size;
    const int numTransitionCounts = 10;
    std::vector<double>& expectedCounts = extras.expectedCounts;
    expectedCounts.resize( numEmissionCounts + numTransitionCounts );
    countSeedMatches( &expectedCounts[0],
		      seq1 + seed.beg1(), seq1 + seed.end1(),
		      seq2 + seed.beg2(), alph );
    expectedCounts[ numEmissionCounts ] += seed.size;  // match count
  }

  // extend a gapped alignment in the left/reverse direction from the seed:
  std::vector<char>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  extend( blocks, columnAmbiguityCodes, centroid, greedyAligner, isGreedy,
	  seq1, seq2, seed.beg1(), seed.beg2(), false, globality,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // maybe unnecessary?

  // convert left-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedBeg1 = seed.beg1();
  SegmentPair::indexT seedBeg2 = aaToDna( seed.beg2(), frameSize );
  for( IT(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    i->start1 = seedBeg1 - i->start1 - i->size;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedBeg2 - i->start2, frameSize ) - i->size;
  }

  // extend a gapped alignment in the right/forward direction from the seed:
  std::vector<SegmentPair> forwardBlocks;
  std::vector<char> forwardAmbiguities;
  extend( forwardBlocks, forwardAmbiguities, centroid, greedyAligner, isGreedy,
	  seq1, seq2, seed.end1(), seed.end2(), true, globality,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // maybe unnecessary?

  // convert right-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedEnd1 = seed.end1();
  SegmentPair::indexT seedEnd2 = aaToDna( seed.end2(), frameSize );
  for( IT(SegmentPair) i = forwardBlocks.begin(); i < forwardBlocks.end();
       ++i ){
    i->start1 = seedEnd1 + i->start1;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedEnd2 + i->start2, frameSize );
  }

  bool isMergeSeedReverse = !blocks.empty() && isNext( blocks.back(), seed );
  bool isMergeSeedForward =
    !forwardBlocks.empty() && isNext( seed, forwardBlocks.back() );

  if( seed.size == 0 && !isMergeSeedReverse && !isMergeSeedForward ){
    // unusual, weird case: give up
    score = -INF;
    return;
  }

  // splice together the two extensions and the seed (a bit messy):

  blocks.reserve( blocks.size() + forwardBlocks.size() +
		  1 - isMergeSeedReverse - isMergeSeedForward );

  if( isMergeSeedReverse ) blocks.back().size += seed.size;
  else                     blocks.push_back(seed);

  if( isMergeSeedForward ){
    blocks.back().size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );

  if( outputType > 3 ){  // set the un-ambiguity of the core to a max value:
    columnAmbiguityCodes.insert( columnAmbiguityCodes.end(), seed.size, 126 );
  }

  columnAmbiguityCodes.insert( columnAmbiguityCodes.end(),
                               forwardAmbiguities.rbegin(),
                               forwardAmbiguities.rend() );
}

// cost of the gap between x and y
static int gapCost(const SegmentPair &x, const SegmentPair &y,
		   const mcf::GapCosts &gapCosts,
		   int frameshiftCost, size_t frameSize) {
  size_t gapSize1 = y.beg1() - x.end1();
  size_t gapSize2, frameshift2;
  sizeAndFrameshift(x.end2(), y.beg2(), frameSize, gapSize2, frameshift2);
  int cost = gapCosts.cost(gapSize1, gapSize2);
  if (frameshift2) cost += frameshiftCost;
  return cost;
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix, int maxDrop,
			   const mcf::GapCosts& gapCosts,
			   int frameshiftCost, size_t frameSize,
			   const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2 ) const{
  int maxScore = 0;
  int runningScore = 0;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];

    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];
      runningScore -= gapCost( x, y, gapCosts, frameshiftCost, frameSize );
      if( !globality && runningScore <= 0 ) return false;
      if( runningScore < maxScore - maxDrop ) return false;
    }

    const uchar* s1 = seq1 + y.beg1();
    const uchar* s2 = seq2 + y.beg2();
    const uchar* e1 = seq1 + y.end1();

    const ScoreMatrixRow* p2 = pssm2 ? pssm2 + y.beg2() : 0;
    const uchar* q1 = qual1 ? qual1 + y.beg1() : 0;
    const uchar* q2 = qual2 ? qual2 + y.beg2() : 0;

    while( s1 < e1 ){
      /**/ if( sm2qual ) runningScore += sm2qual( *s1++, *s2++, *q1++, *q2++ );
      else if( pssm2 )   runningScore += ( *p2++ )[ *s1++ ];
      else               runningScore += scoreMatrix[ *s1++ ][ *s2++ ];

      if( runningScore > maxScore ) maxScore = runningScore;
      else if( !globality && runningScore <= 0 ) return false;
      else if( !globality && (s1 == e1 && i+1 == blocks.size()) ) return false;
      else if( runningScore < maxScore - maxDrop ) return false;
    }
  }

  return true;
}

bool Alignment::hasGoodSegment(const uchar *seq1, const uchar *seq2,
			       int minScore, const ScoreMatrixRow *scoreMatrix,
			       const mcf::GapCosts &gapCosts,
			       int frameshiftCost, size_t frameSize,
			       const ScoreMatrixRow *pssm2,
			       const TwoQualityScoreMatrix &sm2qual,
			       const uchar *qual1, const uchar *qual2) const {
  int score = 0;

  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair& y = blocks[i];

    if (i > 0) {  // between each pair of aligned blocks:
      score -= gapCost(blocks[i - 1], y, gapCosts, frameshiftCost, frameSize);
      if (score < 0) score = 0;
    }

    const uchar *s1 = seq1 + y.beg1();
    const uchar *s2 = seq2 + y.beg2();
    const uchar *e1 = seq1 + y.end1();

    const ScoreMatrixRow *p2 = pssm2 ? pssm2 + y.beg2() : 0;
    const uchar *q1 = qual1 ? qual1 + y.beg1() : 0;
    const uchar *q2 = qual2 ? qual2 + y.beg2() : 0;

    while (s1 < e1) {
      /**/ if (sm2qual) score += sm2qual(*s1++, *s2++, *q1++, *q2++);
      else if (pssm2)   score += (*p2++)[*s1++];
      else              score += scoreMatrix[*s1++][*s2++];

      if (score >= minScore) return true;
      if (score < 0) score = 0;
    }
  }

  return false;
}

static void getColumnCodes(const Centroid& centroid, std::vector<char>& codes,
			   const std::vector<SegmentPair>& chunks,
			   bool isForward) {
  for (size_t i = 0; i < chunks.size(); ++i) {
    const SegmentPair& x = chunks[i];
    centroid.getMatchAmbiguities(codes, x.end1(), x.end2(), x.size);
    size_t j = i + 1;
    bool isNext = (j < chunks.size());
    size_t end1 = isNext ? chunks[j].end1() : 0;
    size_t end2 = isNext ? chunks[j].end2() : 0;
    // ASSUMPTION: if there is an insertion adjacent to a deletion,
    // the deletion will get printed first.
    if (isForward) {
      centroid.getInsertAmbiguities(codes, x.beg2(), end2);
      centroid.getDeleteAmbiguities(codes, x.beg1(), end1);
    } else {
      centroid.getDeleteAmbiguities(codes, x.beg1(), end1);
      centroid.getInsertAmbiguities(codes, x.beg2(), end2);
    }
  }
}

void Alignment::extend( std::vector< SegmentPair >& chunks,
			std::vector< char >& columnCodes,
			Centroid& centroid,
			GreedyXdropAligner& greedyAligner, bool isGreedy,
			const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2,
			bool isForward, int globality,
			const ScoreMatrixRow* sm, int smMax, int maxDrop,
			const mcf::GapCosts& gap,
			int frameshiftCost, size_t frameSize,
			const ScoreMatrixRow* pssm2,
                        const TwoQualityScoreMatrix& sm2qual,
                        const uchar* qual1, const uchar* qual2,
			const Alphabet& alph, AlignmentExtras& extras,
			double gamma, int outputType ){
  const mcf::GapCosts::Piece &del = gap.delPieces[0];
  const mcf::GapCosts::Piece &ins = gap.insPieces[0];
  GappedXdropAligner& aligner = centroid.aligner();

  if( frameSize ){
    assert( outputType < 4 );
    assert( !isGreedy );
    assert( !globality );
    assert( !pssm2 );
    assert( !sm2qual );

    size_t dnaStart = aaToDna( start2, frameSize );
    size_t f = dnaStart + 1;
    size_t r = dnaStart - 1;
    size_t frame1 = dnaToAa( isForward ? f : r, frameSize );
    size_t frame2 = dnaToAa( isForward ? r : f, frameSize );

    score += aligner.align3( seq1 + start1, seq2 + start2,
			     seq2 + frame1, seq2 + frame2, isForward,
			     sm, del.openCost, del.growCost, gap.pairCost,
			     frameshiftCost, maxDrop, smMax );

    size_t end1, end2, size;
    // This should be OK even if end2 < size * 3:
    while( aligner.getNextChunk3( end1, end2, size,
				  del.openCost, del.growCost, gap.pairCost,
				  frameshiftCost ) )
      chunks.push_back( SegmentPair( end1 - size, end2 - size * 3, size ) );

    return;
  }

  int extensionScore =
    isGreedy  ? greedyAligner.align( seq1 + start1, seq2 + start2,
				     isForward, sm, maxDrop, alph.size )
    : sm2qual ? aligner.align2qual( seq1 + start1, qual1 + start1,
				    seq2 + start2, qual2 + start2,
				    isForward, globality, sm2qual,
				    del.openCost, del.growCost,
				    ins.openCost, ins.growCost,
				    gap.pairCost, maxDrop, smMax )
    : pssm2   ? aligner.alignPssm( seq1 + start1, pssm2 + start2,
				   isForward, globality,
				   del.openCost, del.growCost,
				   ins.openCost, ins.growCost,
				   gap.pairCost, maxDrop, smMax )
    :           aligner.align( seq1 + start1, seq2 + start2,
			       isForward, globality, sm,
			       del.openCost, del.growCost,
			       ins.openCost, ins.growCost,
			       gap.pairCost, maxDrop, smMax );

  if( extensionScore == -INF ){
    score = -INF;  // avoid score overflow
    return;  // avoid ill-defined probabilistic alignment
  }

  score += extensionScore;

  if( outputType < 5 || outputType > 6 ){  // ordinary max-score alignment
    size_t end1, end2, size;
    if( isGreedy ){
      while( greedyAligner.getNextChunk( end1, end2, size ) )
	chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
    }else{
      while( aligner.getNextChunk( end1, end2, size,
				   del.openCost, del.growCost,
				   ins.openCost, ins.growCost, gap.pairCost ) )
	chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
    }
  }

  if( outputType > 3 ){  // calculate match probabilities
    assert( !isGreedy );
    assert( !sm2qual );
    if (!isForward) {
      --start1;
      --start2;
    }
    centroid.doForwardBackwardAlgorithm(seq1, seq2, start1, start2, isForward,
					gap, globality);

    if( outputType > 4 && outputType < 7 ){  // gamma-centroid / LAMA alignment
      centroid.dp( gamma );
      centroid.traceback( chunks, gamma );
    }

    getColumnCodes(centroid, columnCodes, chunks, isForward);
    extras.fullScore += centroid.logPartitionFunction();

    if( outputType == 7 ){
      ExpectedCount ec;
      centroid.computeExpectedCounts( seq1, seq2, start1, start2,
				      isForward, gap, alph.size, ec );
      addExpectedCounts( &extras.expectedCounts[0], ec, alph );
    }
  }
}
