// Copyright 2008, 2009, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#ifndef CENTROID_HH
#define CENTROID_HH
#include "GappedXdropAligner.hh"
#include "mcf_gap_costs.hh"
#include "SegmentPair.hh"
#include "OneQualityScoreMatrix.hh"
#include <stddef.h>  // size_t
#include <vector>
#include <iostream> // for debug

namespace cbrc{

  struct ExpectedCount{
  public:
    double emit[scoreMatrixRowSize][scoreMatrixRowSize];
    double toMatch;
    double MD, MP, MI;
    double DD, DI;
    double PP, PD, PI;
    double II;
  public:
    ExpectedCount ();
  };

  /**
   * (1) Forward and backward algorithm on the DP region given by Xdrop algorithm
   * (2) \gamma-centroid decoding
   */
  class Centroid{
  public:
    GappedXdropAligner& aligner() { return xa; }

    // Setters
    void setScoreMatrix( const ScoreMatrixRow* sm, double T );
    void setPssm ( const ScoreMatrixRow* pssm, size_t qsize, double T,
                   const OneQualityExpMatrix& oqem,
                   const uchar* sequenceBeg, const uchar* qualityBeg );
    void setOutputType( int m ) { outputType = m; }

    // For a sequence with quality data, store the probability that
    // each position is each letter (possibly scaled by a constant per
    // position).  xxx I don't think this really belongs in Centroid.
    void setLetterProbsPerPosition(unsigned alphabetSize,
				   size_t sequenceLength,
				   const uchar *sequence,
				   const uchar *qualityCodes,
				   bool isFastq,
				   const double *qualToProbCorrect,
				   const double *letterProbs,
				   const uchar *toUnmasked);

    // start1 is the index of the first letter to look at in seq1
    // start2 is the index of the first letter to look at in seq2

    void doForwardBackwardAlgorithm(const uchar* seq1, const uchar* seq2,
				    size_t start1, size_t start2,
				    bool isExtendFwd, const mcf::GapCosts& gap,
				    int globality) {
      seq1 += start1;
      seq2 += start2;
      const ExpMatrixRow *pssm = isPssm ? pssmExp2 + start2 : 0;
      numAntidiagonals = xa.numAntidiagonals();
      scale.assign(numAntidiagonals + 2, 1.0);
      forward(seq1, seq2, pssm, isExtendFwd, gap, globality);
      mD.assign(numAntidiagonals + 2, 0.0);
      mI.assign(numAntidiagonals + 2, 0.0);
      backward(seq1, seq2, pssm, isExtendFwd, gap, globality);
    }

    double dp( double gamma );

    void traceback(std::vector<SegmentPair> &chunks, double gamma) const {
      if (outputType==5) traceback_centroid(chunks, gamma);
      if (outputType==6) traceback_ama(chunks, gamma);
    }

    double dp_centroid( double gamma );
    void traceback_centroid( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_ama( double gamma );
    void traceback_ama( std::vector< SegmentPair >& chunks, double gamma ) const;

    void getMatchAmbiguities(std::vector<char>& ambiguityCodes,
			     size_t seq1end, size_t seq2end,
			     size_t length) const;

    void getDeleteAmbiguities(std::vector<char>& ambiguityCodes,
			      size_t seq1end, size_t seq1beg) const;

    void getInsertAmbiguities(std::vector<char>& ambiguityCodes,
			      size_t seq2end, size_t seq2beg) const;

    double logPartitionFunction() const;  // a.k.a. full score, forward score

    // Added by MH (2008/10/10) : compute expected counts for transitions and emissions
    void computeExpectedCounts(const uchar* seq1, const uchar* seq2,
			       size_t start1, size_t start2, bool isExtendFwd,
			       const mcf::GapCosts& gap, unsigned alphabetSize,
			       ExpectedCount& count) const;

  private:
    typedef double ExpMatrixRow[scoreMatrixRowSize];

    GappedXdropAligner xa;
    double T; // temperature
    size_t numAntidiagonals;
    double match_score[scoreMatrixRowSize][scoreMatrixRowSize];
    bool isPssm;
    std::vector<double> pssmExp; //
    ExpMatrixRow* pssmExp2; // pre-computed pssm for prob align
    std::vector<double> letterProbsPerPosition;  // for uncertain sequences
    int outputType;

    typedef std::vector< double > dvec_t;

    dvec_t fM; // f^M(i,j)
    dvec_t fD; // f^D(i,j) Ix
    dvec_t fI; // f^I(i,j) Iy
    dvec_t fP; // f^P(i,j)

    dvec_t bM; // b^M(i,j)
    dvec_t bD; // b^D(i,j)
    dvec_t bI; // b^I(i,j)
    dvec_t bP; // b^P(i,j)

    dvec_t mD;
    dvec_t mI;
    dvec_t mX1;
    dvec_t mX2;

    dvec_t X; // DP tables for $gamma$-decoding

    dvec_t scale; // scale[n] is a scaling factor for the n-th anti-diagonal

    double bestScore;
    size_t bestAntiDiagonal;
    size_t bestPos1;

    void forward(const uchar* seq1, const uchar* seq2,
		 const ExpMatrixRow* pssm, bool isExtendFwd,
		 const mcf::GapCosts& gap, int globality);

    void backward(const uchar* seq1, const uchar* seq2,
		  const ExpMatrixRow* pssm, bool isExtendFwd,
		  const mcf::GapCosts& gap, int globality);

    void initForwardMatrix();
    void initBackwardMatrix();

    void updateScore(double score, size_t antiDiagonal, size_t cur) {
      if (bestScore < score) {
	bestScore = score;
	bestAntiDiagonal = antiDiagonal;
	bestPos1 = cur;
      }
    }

    // start of the x-drop region (i.e. number of skipped seq1 letters
    // before the x-drop region) for this antidiagonal
    size_t seq1start( size_t antidiagonal ) const {
      return xa.scoreEndIndex( antidiagonal ) - xa.scoreOrigin( antidiagonal );
    }
  };

}

#endif
