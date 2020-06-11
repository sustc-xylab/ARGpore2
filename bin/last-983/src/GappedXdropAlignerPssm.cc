// Copyright 2011, 2012, 2013 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

namespace cbrc {

int GappedXdropAligner::alignPssm(const uchar *seq,
                                  const ScoreMatrixRow *pssm,
                                  bool isForward,
				  int globality,
				  int delExistenceCost,
				  int delExtensionCost,
				  int insExistenceCost,
				  int insExtensionCost,
                                  int gapUnalignedCost,
                                  int maxScoreDrop,
                                  int maxMatchScore) {
  const bool isAffine = isAffineGaps(delExistenceCost, delExtensionCost,
				     insExistenceCost, insExtensionCost,
				     gapUnalignedCost);

  std::size_t maxSeq1begs[] = { 0, 9 };
  std::size_t minSeq1ends[] = { 1, 0 };

  int bestScore = 0;
  int bestEdgeScore = -INF;
  std::size_t bestEdgeAntidiagonal = 0;
  std::size_t bestEdgeSeq1position = 0;

  init();

  for (std::size_t antidiagonal = 0; /* noop */; ++antidiagonal) {
    std::size_t seq1beg = std::min(maxSeq1begs[0], maxSeq1begs[1]);
    std::size_t seq1end = std::max(minSeq1ends[0], minSeq1ends[1]);

    if (seq1beg >= seq1end) break;

    std::size_t scoreEnd = scoreEnds.back();
    std::size_t numCells = seq1end - seq1beg;

    initAntidiagonal(seq1beg, scoreEnd, numCells);

    std::size_t seq2pos = antidiagonal - seq1beg;

    const uchar *s1 = isForward ? seq + seq1beg : seq - seq1beg - 1;
    const ScoreMatrixRow *s2 = isForward ? pssm + seq2pos : pssm - seq2pos - 1;

    if (!globality && isDelimiter(0, *s2))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y1 = &yScores[hori(antidiagonal, seq1beg)];
    const int *z1 = &zScores[vert(antidiagonal, seq1beg)];
    const int *x2 = &xScores[diag(antidiagonal, seq1beg)];

    const int *x0last = x0 + numCells;

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0base = x0 - seq1beg;

    if (globality && isDelimiter(0, *s2)) {
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      int b = maxValue(*x2, *z1 - insExtensionCost, *z2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestEdgeSeq1position,
		    b, antidiagonal, seq1beg);
    }

    if (isAffine) {
      if (isForward)
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + (*s2)[*s1];
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  --s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
      else
        while (1) {
          int x = *x2;
          int y = *y1 - delExtensionCost;
          int z = *z1 - insExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + (*s2)[*s1];
            *y0 = maxValue(b - delExistenceCost, y);
            *z0 = maxValue(b - insExistenceCost, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  ++s2;  ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;
        }
    } else {
      const int *y2 = &yScores[diag(antidiagonal, seq1beg)];
      const int *z2 = &zScores[diag(antidiagonal, seq1beg)];
      while (1) {
        int x = *x2;
        int y = maxValue(*y1 - delExtensionCost, *y2 - gapUnalignedCost);
        int z = maxValue(*z1 - insExtensionCost, *z2 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
          updateBest(bestScore, b, antidiagonal, x0, x0base);
          *x0 = b + (*s2)[*s1];
          *y0 = maxValue(b - delExistenceCost, y);
          *z0 = maxValue(b - insExistenceCost, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y1;  ++z1;  ++x2;  ++y2;  ++z2;
        if (isForward) { ++s1;  --s2; }
        else           { --s1;  ++s2; }
      }
    }

    if (globality && isDelimiter(*s1, *pssm)) {
      const int *y2 = &yScores[diag(antidiagonal, seq1end-1)];
      int b = maxValue(*x2, *y1 - delExtensionCost, *y2 - gapUnalignedCost);
      if (b >= minScore)
	updateBest1(bestEdgeScore, bestEdgeAntidiagonal, bestEdgeSeq1position,
		    b, antidiagonal, seq1end-1);
    }

    if (!globality && isDelimiter(*s1, *pssm))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    updateFiniteEdges(maxSeq1begs, minSeq1ends, x0base, x0 + 1, numCells);
  }

  if (globality) {
    bestAntidiagonal = bestEdgeAntidiagonal;
    bestSeq1position = bestEdgeSeq1position;
    bestScore = bestEdgeScore;
  }
  return bestScore;
}

}
