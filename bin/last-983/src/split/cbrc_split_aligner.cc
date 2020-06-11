// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

#include "cbrc_split_aligner.hh"
#include "mcf_substitution_matrix_stats.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstring>  // strcmp
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

namespace cbrc {

template<typename T, int N> T arrayMax(T (&array)[N]) {
  return *std::max_element(array, array + N);
}

// Orders candidate alignments by increasing DP start coordinate.
// Breaks ties by chromosome & strand, then by increasing genomic
// start coordinate.
struct QbegLess {
  QbegLess(const unsigned *b, const unsigned *r, const unsigned *rb)
    : dpBegs(b), rnameAndStrandIds(r), rBegs(rb) {}

  bool operator()(unsigned a, unsigned b) const {
    return dpBegs[a] != dpBegs[b]
      ? dpBegs[a] < dpBegs[b]
      : rnameAndStrandIds[a] != rnameAndStrandIds[b]
      ? rnameAndStrandIds[a] < rnameAndStrandIds[b]
      : rBegs[a] < rBegs[b];
  }

  const unsigned *dpBegs;
  const unsigned *rnameAndStrandIds;
  const unsigned *rBegs;
};

// Orders candidate alignments by decreasing DP end coordinate.
// Breaks ties by chromosome & strand, then by decreasing genomic end
// coordinate.
struct QendLess {
  QendLess(const unsigned *e, const unsigned *r, const unsigned *re)
    : dpEnds(e), rnameAndStrandIds(r), rEnds(re) {}

  bool operator()(unsigned a, unsigned b) const {
    return dpEnds[a] != dpEnds[b]
      ? dpEnds[a] > dpEnds[b]
      : rnameAndStrandIds[a] != rnameAndStrandIds[b]
      ? rnameAndStrandIds[a] < rnameAndStrandIds[b]
      : rEnds[a] > rEnds[b];
  }

  const unsigned *dpEnds;
  const unsigned *rnameAndStrandIds;
  const unsigned *rEnds;
};

// Orders candidate alignments by: chromosome & strand, then
// increasing genomic start coordinate.
struct RbegLess {
  RbegLess(const unsigned *r, const unsigned *rb)
    : rnameAndStrandIds(r), rBegs(rb) {}

  bool operator()(unsigned a, unsigned b) const {
    return rnameAndStrandIds[a] != rnameAndStrandIds[b]
      ? rnameAndStrandIds[a] < rnameAndStrandIds[b]
      : rBegs[a] < rBegs[b];
  }

  const unsigned *rnameAndStrandIds;
  const unsigned *rBegs;
};

// Orders candidate alignments by: chromosome & strand, then
// decreasing genomic end coordinate.
struct RendLess {
  RendLess(const unsigned *r, const unsigned *re)
    : rnameAndStrandIds(r), rEnds(re) {}

  bool operator()(unsigned a, unsigned b) const {
    return rnameAndStrandIds[a] != rnameAndStrandIds[b]
      ? rnameAndStrandIds[a] < rnameAndStrandIds[b]
      : rEnds[a] > rEnds[b];
  }

  const unsigned *rnameAndStrandIds;
  const unsigned *rEnds;
};

// Merges the elements of the sorted range [beg2,end2) into the sorted
// range [beg1,end1).  Assumes it's OK to add elements past end1.
template<typename T>
void mergeInto(unsigned* beg1,
	       unsigned* end1,
	       const unsigned* beg2,
	       const unsigned* end2,
	       T lessFunc) {
  unsigned* end3 = end1 + (end2 - beg2);
  for (;;) {
    if (beg2 == end2)
      break;
    if (beg1 == end1) {
      std::copy(beg2, end2, beg1);
      break;
    }
    --end3;
    if (lessFunc(*(end2-1), *(end1-1)))
      *end3 = *--end1;
    else
      *end3 = *--end2;
  }
}

int SplitAligner::score_mat[64][64][numQualCodes];

// The score for a cis-splice with the given distance (i.e. intron length)
int SplitAligner::calcSpliceScore(double dist) const {
    double logDist = std::log(dist);
    double d = logDist - meanLogDist;
    double s = spliceTerm1 + spliceTerm2 * d * d - logDist;
    return std::floor(scale * s + 0.5);
}

// The dinucleotide immediately downstream on the forward strand
static unsigned spliceBegSignalFwd(const uchar *seqPtr,
				   const uchar *toUnmasked) {
  unsigned n1 = toUnmasked[*seqPtr];
  if (n1 >= 4) return 16;
  unsigned n2 = toUnmasked[*(seqPtr + 1)];
  if (n2 >= 4) return 16;
  return n1 * 4 + n2;
}

// The dinucleotide immediately downstream on the reverse strand
static unsigned spliceBegSignalRev(const uchar *seqPtr,
				   const uchar *toUnmasked) {
  unsigned n1 = toUnmasked[*(seqPtr - 1)];
  if (n1 >= 4) return 16;
  unsigned n2 = toUnmasked[*(seqPtr - 2)];
  if (n2 >= 4) return 16;
  return 15 - (n1 * 4 + n2);  // reverse-complement
}

// The dinucleotide immediately upstream on the forward strand
static unsigned spliceEndSignalFwd(const uchar *seqPtr,
				   const uchar *toUnmasked) {
  unsigned n2 = toUnmasked[*(seqPtr - 1)];
  if (n2 >= 4) return 16;
  unsigned n1 = toUnmasked[*(seqPtr - 2)];
  if (n1 >= 4) return 16;
  return n1 * 4 + n2;
}

// The dinucleotide immediately upstream on the reverse strand
static unsigned spliceEndSignalRev(const uchar *seqPtr,
				   const uchar *toUnmasked) {
  unsigned n2 = toUnmasked[*seqPtr];
  if (n2 >= 4) return 16;
  unsigned n1 = toUnmasked[*(seqPtr + 1)];
  if (n1 >= 4) return 16;
  return 15 - (n1 * 4 + n2);  // reverse-complement
}

unsigned SplitAligner::findScore(unsigned j, long score) const {
  for (unsigned i = 0; i < numAlns; ++i) {
    if (dpBeg(i) >= j || dpEnd(i) < j) continue;
    size_t ij = matrixRowOrigins[i] + j;
    if (Vmat[ij] + spliceBegScore(ij) == score) return i;
  }
  return numAlns;
}

unsigned SplitAligner::findSpliceScore(unsigned i, unsigned j,
				       long score) const {
    assert(splicePrior > 0.0);
    size_t ij = matrixRowOrigins[i] + j;
    unsigned iSeq = rnameAndStrandIds[i];
    unsigned iEnd = spliceEndCoords[ij];
    int iScore = spliceEndScore(ij);
    for (unsigned k = 0; k < numAlns; k++) {
	if (rnameAndStrandIds[k] != iSeq) continue;
	if (dpBeg(k) >= j || dpEnd(k) < j) continue;
	size_t kj = matrixRowOrigins[k] + j;
	unsigned kBeg = spliceBegCoords[kj];
	if (iEnd <= kBeg) continue;
	int s = iScore + spliceBegScore(kj) + spliceScore(iEnd - kBeg);
	if (Vmat[kj] + s == score) return k;
    }
    return numAlns;
}

long SplitAligner::scoreFromSplice(unsigned i, unsigned j,
				   unsigned oldNumInplay,
				   unsigned& oldInplayPos) const {
  size_t ij = matrixRowOrigins[i] + j;
  long score = LONG_MIN;
  unsigned iSeq = rnameAndStrandIds[i];
  unsigned iEnd = spliceEndCoords[ij];

  for (/* noop */; oldInplayPos < oldNumInplay; ++oldInplayPos) {
    unsigned k = oldInplayAlnIndices[oldInplayPos];
    if (rnameAndStrandIds[k] < iSeq) continue;
    if (rnameAndStrandIds[k] > iSeq || rBegs[k] >= iEnd) return score;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kBeg = spliceBegCoords[kj];
    if (kBeg >= rBegs[i] || rBegs[i] - kBeg <= maxSpliceDist) break;
  }

  for (unsigned y = oldInplayPos; y < oldNumInplay; ++y) {
    unsigned k = oldInplayAlnIndices[y];
    if (rnameAndStrandIds[k] > iSeq || rBegs[k] >= iEnd) break;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kBeg = spliceBegCoords[kj];
    if (iEnd <= kBeg) continue;
    if (iEnd - kBeg > maxSpliceDist) continue;
    score = std::max(score,
		     Vmat[kj] + spliceBegScore(kj) + spliceScore(iEnd - kBeg));
  }

  return score;
}

void SplitAligner::updateInplayAlnIndicesF(unsigned& sortedAlnPos,
					   unsigned& oldNumInplay,
                                           unsigned& newNumInplay,
                                           unsigned j) {  // query coordinate
  oldInplayAlnIndices.swap(newInplayAlnIndices);
  oldNumInplay = newNumInplay;

  unsigned *newBeg = &newInplayAlnIndices[0];
  unsigned *newEnd = newBeg;
  const unsigned *oldBeg = &oldInplayAlnIndices[0];
  const unsigned *oldEnd = oldBeg + oldNumInplay;

  while (oldBeg < oldEnd) {
    unsigned i = *oldBeg++;
    if (dpEnd(i) == j) continue;  // it is no longer "in play"
    *newEnd++ = i;
  }

  unsigned sortedAlnOldPos = sortedAlnPos;

  for (/* noop */; sortedAlnPos < numAlns; ++sortedAlnPos) {
    unsigned i = sortedAlnIndices[sortedAlnPos];
    if (dpBeg(i) > j) break;  // it is not yet "in play"
  }

  mergeInto(newBeg,
            newEnd,
            &sortedAlnIndices[0] + sortedAlnOldPos,
            &sortedAlnIndices[0] + sortedAlnPos,
            RbegLess(&rnameAndStrandIds[0], &rBegs[0]));

  newNumInplay = (newEnd - newBeg) + (sortedAlnPos - sortedAlnOldPos);
}

void SplitAligner::updateInplayAlnIndicesB(unsigned& sortedAlnPos,
					   unsigned& oldNumInplay,
                                           unsigned& newNumInplay,
                                           unsigned j) {  // query coordinate
  oldInplayAlnIndices.swap(newInplayAlnIndices);
  oldNumInplay = newNumInplay;

  unsigned *newBeg = &newInplayAlnIndices[0];
  unsigned *newEnd = newBeg;
  const unsigned *oldBeg = &oldInplayAlnIndices[0];
  const unsigned *oldEnd = oldBeg + oldNumInplay;

  while (oldBeg < oldEnd) {
    unsigned i = *oldBeg++;
    if (dpBeg(i) == j) continue;  // it is no longer "in play"
    *newEnd++ = i;
  }

  unsigned sortedAlnOldPos = sortedAlnPos;

  for (/* noop */; sortedAlnPos < numAlns; ++sortedAlnPos) {
    unsigned i = sortedAlnIndices[sortedAlnPos];
    if (dpEnd(i) < j) break;  // it is not yet "in play"
  }

  mergeInto(newBeg,
            newEnd,
            &sortedAlnIndices[0] + sortedAlnOldPos,
            &sortedAlnIndices[0] + sortedAlnPos,
            RendLess(&rnameAndStrandIds[0], &rEnds[0]));

  newNumInplay = (newEnd - newBeg) + (sortedAlnPos - sortedAlnOldPos);
}

long SplitAligner::viterbi() {
    resizeMatrix(Vmat);
    resizeVector(Vvec);

    for (unsigned i = 0; i < numAlns; ++i) cell(Vmat, i, dpBeg(i)) = INT_MIN/2;
    long maxScore = 0;
    cell(Vvec, minBeg) = maxScore;
    long scoreFromJump = restartScore;

    stable_sort(sortedAlnIndices.begin(), sortedAlnIndices.end(),
		QbegLess(&dpBegs[0], &rnameAndStrandIds[0], &rBegs[0]));
    unsigned sortedAlnPos = 0;
    unsigned oldNumInplay = 0;
    unsigned newNumInplay = 0;

    for (unsigned j = minBeg; j < maxEnd; j++) {
	updateInplayAlnIndicesF(sortedAlnPos, oldNumInplay, newNumInplay, j);
	unsigned oldInplayPos = 0;
	long sMax = INT_MIN/2;
	for (unsigned x = 0; x < newNumInplay; ++x) {
	    unsigned i = newInplayAlnIndices[x];
	    size_t ij = matrixRowOrigins[i] + j;

	    long s = scoreFromJump;
	    if (splicePrior > 0.0)
	      s = std::max(s,
			   scoreFromSplice(i, j, oldNumInplay, oldInplayPos));
	    s += spliceEndScore(ij);
	    s = std::max(s, Vmat[ij] + Dmat[ij]);
	    if (restartProb <= 0 && alns[i].qstart == j && s < 0) s = 0;
	    s += Amat[ij];

	    Vmat[ij + 1] = s;
	    sMax = std::max(sMax, s + spliceBegScore(ij + 1));
	}
	maxScore = std::max(sMax, maxScore);
	cell(Vvec, j+1) = maxScore;
	scoreFromJump = std::max(sMax + jumpScore, maxScore + restartScore);
    }

    return (restartProb > 0) ? maxScore : endScore();
}

long SplitAligner::endScore() const {
    long score = LONG_MIN;
    for (unsigned i = 0; i < numAlns; ++i)
	score = std::max(score, cell(Vmat, i, alns[i].qend));
    return score;
}

unsigned SplitAligner::findEndScore(long score) const {
    for (unsigned i = 0; i < numAlns; ++i)
        if (cell(Vmat, i, alns[i].qend) == score)
            return i;
    return numAlns;
}

void SplitAligner::traceBack(long viterbiScore,
			     std::vector<unsigned>& alnNums,
			     std::vector<unsigned>& queryBegs,
			     std::vector<unsigned>& queryEnds) const {
  unsigned i, j;
  if (restartProb > 0) {
    j = maxEnd;
    long t = cell(Vvec, j);
    if (t == 0) return;
    while (t == cell(Vvec, j-1)) --j;
    i = findScore(j, t);
    assert(i < numAlns);
  } else {
    i = findEndScore(viterbiScore);
    assert(i < numAlns);
    j = alns[i].qend;
  }

  alnNums.push_back(i);
  queryEnds.push_back(j);

  for (;;) {
    --j;
    size_t ij = matrixRowOrigins[i] + j;
    long score = Vmat[ij + 1] - Amat[ij];
    if (restartProb <= 0 && alns[i].qstart == j && score == 0) {
      queryBegs.push_back(j);
      return;
    }

    // We either stay in this alignment, or jump to another one.  If
    // the scores are equally good, then we stay if the strand is "+",
    // else jump.  This gives cleaner inversion boundaries, but it
    // makes some other kinds of boundary less clean.  What's the best
    // procedure for tied scores?

    bool isStay = (score == Vmat[ij] + Dmat[ij]);
    if (isStay && alns[i].isForwardStrand()) continue;

    long s = score - spliceEndScore(ij);
    long t = s - restartScore;
    if (t == cell(Vvec, j)) {
      queryBegs.push_back(j);
      if (t == 0) return;
      while (t == cell(Vvec, j-1)) --j;
      i = findScore(j, t);
    } else {
      if (isStay) continue;
      queryBegs.push_back(j);
      unsigned k = findScore(j, s - jumpScore);
      i = (k < numAlns) ? k : findSpliceScore(i, j, score);
    }
    assert(i < numAlns);
    alnNums.push_back(i);
    queryEnds.push_back(j);
  }
}

int SplitAligner::segmentScore(unsigned alnNum,
			       unsigned queryBeg, unsigned queryEnd) const {
  int score = 0;
  unsigned i = alnNum;
  for (unsigned j = queryBeg; j < queryEnd; ++j) {
    size_t ij = matrixRowOrigins[i] + j;
    score += Amat[ij];
    if (j > queryBeg) score += Dmat[ij];
  }
  return score;
}

double SplitAligner::probFromSpliceF(unsigned i, unsigned j,
				     unsigned oldNumInplay,
				     unsigned& oldInplayPos) const {
  size_t ij = matrixRowOrigins[i] + j;
  double sum = 0.0;
  unsigned iSeq = rnameAndStrandIds[i];
  unsigned iEnd = spliceEndCoords[ij];

  for (/* noop */; oldInplayPos < oldNumInplay; ++oldInplayPos) {
    unsigned k = oldInplayAlnIndices[oldInplayPos];
    if (rnameAndStrandIds[k] < iSeq) continue;
    if (rnameAndStrandIds[k] > iSeq || rBegs[k] >= iEnd) return sum;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kBeg = spliceBegCoords[kj];
    if (kBeg >= rBegs[i] || rBegs[i] - kBeg <= maxSpliceDist) break;
  }

  for (unsigned y = oldInplayPos; y < oldNumInplay; ++y) {
    unsigned k = oldInplayAlnIndices[y];
    if (rnameAndStrandIds[k] > iSeq || rBegs[k] >= iEnd) break;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kBeg = spliceBegCoords[kj];
    if (iEnd <= kBeg) continue;
    if (iEnd - kBeg > maxSpliceDist) continue;
    sum += Fmat[kj] * spliceBegProb(kj) * spliceProb(iEnd - kBeg);
  }

  return sum;
}

double SplitAligner::probFromSpliceB(unsigned i, unsigned j,
				     unsigned oldNumInplay,
				     unsigned& oldInplayPos) const {
  size_t ij = matrixRowOrigins[i] + j;
  double sum = 0.0;
  unsigned iSeq = rnameAndStrandIds[i];
  unsigned iBeg = spliceBegCoords[ij];

  for (/* noop */; oldInplayPos < oldNumInplay; ++oldInplayPos) {
    unsigned k = oldInplayAlnIndices[oldInplayPos];
    if (rnameAndStrandIds[k] < iSeq) continue;
    if (rnameAndStrandIds[k] > iSeq || rEnds[k] <= iBeg) return sum;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kEnd = spliceEndCoords[kj];
    if (kEnd <= rEnds[i] || kEnd - rEnds[i] <= maxSpliceDist) break;
  }

  for (unsigned y = oldInplayPos; y < oldNumInplay; ++y) {
    unsigned k = oldInplayAlnIndices[y];
    if (rnameAndStrandIds[k] > iSeq || rEnds[k] <= iBeg) break;
    size_t kj = matrixRowOrigins[k] + j;
    unsigned kEnd = spliceEndCoords[kj];
    if (kEnd <= iBeg) continue;
    if (kEnd - iBeg > maxSpliceDist) continue;
    sum += Bmat[kj] * spliceEndProb(kj) * spliceProb(kEnd - iBeg);
  }

  return sum;
}

void SplitAligner::forward() {
    resizeVector(rescales);
    cell(rescales, minBeg) = 1.0;

    resizeMatrix(Fmat);
    for (unsigned i = 0; i < numAlns; ++i) cell(Fmat, i, dpBeg(i)) = 0.0;
    double sumProb = 1;
    double probFromJump = restartProb;
    double begprob = 1.0;
    double zF = 0.0;  // sum of probabilities from the forward algorithm

    stable_sort(sortedAlnIndices.begin(), sortedAlnIndices.end(),
		QbegLess(&dpBegs[0], &rnameAndStrandIds[0], &rBegs[0]));
    unsigned sortedAlnPos = 0;
    unsigned oldNumInplay = 0;
    unsigned newNumInplay = 0;

    for (unsigned j = minBeg; j < maxEnd; j++) {
	updateInplayAlnIndicesF(sortedAlnPos, oldNumInplay, newNumInplay, j);
	unsigned oldInplayPos = 0;
	double r = cell(rescales, j);
	zF /= r;
	double pSum = 0.0;
	double rNew = 1.0;
	for (unsigned x = 0; x < newNumInplay; ++x) {
	    unsigned i = newInplayAlnIndices[x];
	    size_t ij = matrixRowOrigins[i] + j;

	    double p = probFromJump;
	    if (splicePrior > 0.0)
	      p += probFromSpliceF(i, j, oldNumInplay, oldInplayPos);
	    p *= spliceEndProb(ij);
	    p += Fmat[ij] * Dexp[ij];
	    if (restartProb <= 0 && alns[i].qstart == j) p += begprob;
	    p = p * Aexp[ij] / r;

	    Fmat[ij + 1] = p;
	    if (alns[i].qend == j+1) zF += p;
	    pSum += p * spliceBegProb(ij + 1);
	    rNew += p;
        }
        begprob /= r;
        cell(rescales, j+1) = rNew;
	sumProb = pSum + sumProb / r;
	probFromJump = pSum * jumpProb + sumProb * restartProb;
    }

    if (restartProb > 0) zF = sumProb;
    //zF /= cell(rescales, maxEnd);
    cell(rescales, maxEnd) = zF;  // this causes scaled zF to equal 1
}

void SplitAligner::backward() {
    resizeMatrix(Bmat);
    for (unsigned i = 0; i < numAlns; ++i) cell(Bmat, i, dpEnd(i)) = 0.0;
    double sumProb = (restartProb > 0) ? 1 : 0;
    double probFromJump = sumProb;
    double endprob = 1.0;
    //double zB = 0.0;  // sum of probabilities from the backward algorithm

    stable_sort(sortedAlnIndices.begin(), sortedAlnIndices.end(),
		QendLess(&dpEnds[0], &rnameAndStrandIds[0], &rEnds[0]));
    unsigned sortedAlnPos = 0;
    unsigned oldNumInplay = 0;
    unsigned newNumInplay = 0;

    for (unsigned j = maxEnd; j > minBeg; j--) {
	updateInplayAlnIndicesB(sortedAlnPos, oldNumInplay, newNumInplay, j);
	unsigned oldInplayPos = 0;
	double r = cell(rescales, j);
	//zB /= r;
	double pSum = 0.0;
	for (unsigned x = 0; x < newNumInplay; ++x) {
	    unsigned i = newInplayAlnIndices[x];
	    size_t ij = matrixRowOrigins[i] + j;

	    double p = probFromJump;
	    if (splicePrior > 0.0)
	      p += probFromSpliceB(i, j, oldNumInplay, oldInplayPos);
	    p *= spliceBegProb(ij);
	    p += Bmat[ij] * Dexp[ij];
	    if (restartProb <= 0 && alns[i].qend == j) p += endprob;
	    p = p * Aexp[ij - 1] / r;

	    Bmat[ij - 1] = p;
	    //if (alns[i].qstart == j-1) zB += p;
	    pSum += p * spliceEndProb(ij - 1);
        }
        endprob /= r;
	sumProb = pSum * restartProb + sumProb / r;
	probFromJump = pSum * jumpProb + sumProb;
    }

    //if (restartProb > 0) zB = sumProb;
    //zB /= cell(rescales, minBeg);
}

std::vector<double>
SplitAligner::marginalProbs(unsigned queryBeg, unsigned alnNum,
			    unsigned alnBeg, unsigned alnEnd) const {
    std::vector<double> output;
    unsigned i = alnNum;
    unsigned j = queryBeg;
    for (unsigned pos = alnBeg; pos < alnEnd; ++pos) {
	size_t ij = matrixRowOrigins[i] + j;
        if (alns[i].qalign[pos] == '-') {
            double value = Fmat[ij] * Bmat[ij] * Dexp[ij] / cell(rescales, j);
            output.push_back(value);
        } else {
            double value = Fmat[ij + 1] * Bmat[ij] / Aexp[ij];
            if (value != value) value = 0.0;
            output.push_back(value);
            j++;
        }
    }
    return output;
}

// The next routine represents affine gap scores in a cunning way.
// Amat holds scores at query bases, and at every base that is aligned
// to a gap it gets a score of insExistenceScore + insExtensionScore.
// Dmat holds scores between query bases, and between every pair of
// bases that are both aligned to gaps it gets a score of
// -insExistenceScore.  This produces suitable affine gap scores, even
// if we jump from one alignment to another in the middle of a gap.

void SplitAligner::calcBaseScores(unsigned i) {
  const int firstInsScore = insExistenceScore + insExtensionScore;
  const int tweenInsScore = -insExistenceScore;

  const UnsplitAlignment& a = alns[i];
  const size_t origin = matrixRowOrigins[i];

  int *AmatB = &Amat[origin + dpBeg(i)];
  int *DmatB = &Dmat[origin + dpBeg(i)];
  int *AmatS = &Amat[origin + a.qstart];
  int *AmatE = &Amat[origin + dpEnd(i)];

  int delScore = 0;
  int insCompensationScore = 0;

  // treat any query letters before the alignment as insertions:
  while (AmatB < AmatS) {
    *AmatB++ = firstInsScore;
    *DmatB++ = delScore + insCompensationScore;
    delScore = 0;
    insCompensationScore = tweenInsScore;
  }

  const char *rAlign = a.ralign;
  const char *qAlign = a.qalign;
  const char *qQual = a.qQual;

  while (*qAlign) {
    unsigned char x = *rAlign++;
    unsigned char y = *qAlign++;
    int q = qQual ? (*qQual++ - qualityOffset) : (numQualCodes - 1);
    if (x == '-') {  // gap in reference sequence: insertion
      *AmatB++ = firstInsScore;
      *DmatB++ = delScore + insCompensationScore;
      delScore = 0;
      insCompensationScore = tweenInsScore;
    } else if (y == '-') {  // gap in query sequence: deletion
      if (delScore == 0) delScore = gapExistenceScore;
      delScore += gapExtensionScore;
      insCompensationScore = 0;
    } else {
      assert(q >= 0);
      if (q >= numQualCodes) q = numQualCodes - 1;
      *AmatB++ = score_mat[x % 64][y % 64][q];
      *DmatB++ = delScore;
      delScore = 0;
      insCompensationScore = 0;
    }
    // Amazingly, in ASCII, '.' equals 'n' mod 64.
    // So '.' will get the same scores as 'n'.
  }

  // treat any query letters after the alignment as insertions:
  while (AmatB < AmatE) {
    *AmatB++ = firstInsScore;
    *DmatB++ = delScore + insCompensationScore;
    delScore = 0;
    insCompensationScore = tweenInsScore;
  }

  *DmatB++ = delScore;
}

void SplitAligner::initRbegsAndEnds() {
  for (unsigned i = 0; i < numAlns; ++i) {
    const UnsplitAlignment& a = alns[i];
    rBegs[i] = a.rstart;
    rEnds[i] = a.rend;
  }
}

void SplitAligner::initSpliceCoords(unsigned i) {
  const UnsplitAlignment& a = alns[i];
  unsigned j = dpBeg(i);
  unsigned k = a.rstart;

  cell(spliceBegCoords, i, j) = k;
  while (j < a.qstart) {
    cell(spliceEndCoords, i, j) = k;
    ++j;
    cell(spliceBegCoords, i, j) = k;
  }
  for (unsigned x = 0; a.ralign[x]; ++x) {
    if (a.qalign[x] != '-') cell(spliceEndCoords, i, j) = k;
    if (a.qalign[x] != '-') ++j;
    if (a.ralign[x] != '-') ++k;
    if (a.qalign[x] != '-') cell(spliceBegCoords, i, j) = k;
  }
  while (j < dpEnd(i)) {
    cell(spliceEndCoords, i, j) = k;
    ++j;
    cell(spliceBegCoords, i, j) = k;
  }
  cell(spliceEndCoords, i, j) = k;

  assert(k == a.rend);  // xxx
}

static const uchar *seqBeg(const MultiSequence &m, size_t sequenceIndex) {
  return m.seqReader() + m.seqBeg(sequenceIndex);
}

static const uchar *seqEnd(const MultiSequence &m, size_t sequenceIndex) {
  return m.seqReader() + m.seqEnd(sequenceIndex);
}

void SplitAligner::initSpliceSignals(unsigned i) {
  const uchar *toUnmasked = alphabet.numbersToUppercase;
  const UnsplitAlignment& a = alns[i];

  StringNumMap::const_iterator f = chromosomeIndex.find(a.rname);
  if (f == chromosomeIndex.end())
    err("can't find " + std::string(a.rname) + " in the genome");
  size_t v = f->second % maxGenomeVolumes();
  size_t c = f->second / maxGenomeVolumes();
  const uchar *chromBeg = seqBeg(genome[v], c);
  const uchar *chromEnd = seqEnd(genome[v], c);
  if (a.rend > chromEnd - chromBeg)
    err("alignment beyond the end of " + std::string(a.rname));

  size_t rowBeg = matrixRowOrigins[i] + dpBeg(i);
  const unsigned *begCoords = &spliceBegCoords[rowBeg];
  const unsigned *endCoords = &spliceEndCoords[rowBeg];
  unsigned char *begSignals = &spliceBegSignals[rowBeg];
  unsigned char *endSignals = &spliceEndSignals[rowBeg];
  unsigned dpLen = dpEnd(i) - dpBeg(i);

  if (a.isForwardStrand()) {
    for (unsigned j = 0; j <= dpLen; ++j) {
      begSignals[j] = spliceBegSignalFwd(chromBeg + begCoords[j], toUnmasked);
      endSignals[j] = spliceEndSignalFwd(chromBeg + endCoords[j], toUnmasked);
    }
  } else {
    for (unsigned j = 0; j <= dpLen; ++j) {
      begSignals[j] = spliceBegSignalRev(chromEnd - begCoords[j], toUnmasked);
      endSignals[j] = spliceEndSignalRev(chromEnd - endCoords[j], toUnmasked);
    }
  }
}

const uchar sequenceEndSentinel = 4;

static void getNextSignal(uchar *out, const uchar *seq) {
  out[0] = seq[0];
  out[1] = (seq[0] == sequenceEndSentinel) ? sequenceEndSentinel : seq[1];
}

static void getPrevSignal(uchar *out, const uchar *seq) {
  out[1] = seq[-1];
  out[0] = (seq[-1] == sequenceEndSentinel) ? sequenceEndSentinel : seq[-2];
}

static char decodeOneBase(const uchar *decode, uchar x) {
  return (x == sequenceEndSentinel) ? 'N' : decode[x];
}

static void decodeSpliceSignal(char *out,
			       const uchar *signal,
			       const uchar *decode,
			       const uchar *complement,
			       bool isSameStrand) {
  if (isSameStrand) {
    out[0] = decodeOneBase(decode, signal[0]);
    out[1] = decodeOneBase(decode, signal[1]);
  } else {
    out[0] = decodeOneBase(decode, complement[signal[1]]);
    out[1] = decodeOneBase(decode, complement[signal[0]]);
  }
}

void SplitAligner::spliceBegSignal(char *out,
				   unsigned alnNum, unsigned queryPos,
				   bool isSenseStrand) const {
  const UnsplitAlignment& a = alns[alnNum];
  bool isForwardStrand = a.isForwardStrand();
  StringNumMap::const_iterator f = chromosomeIndex.find(a.rname);
  size_t v = f->second % maxGenomeVolumes();
  size_t c = f->second / maxGenomeVolumes();
  uchar signal[2];
  unsigned coord = cell(spliceBegCoords, alnNum, queryPos);
  if (isForwardStrand) getNextSignal(signal, seqBeg(genome[v], c) + coord);
  else                 getPrevSignal(signal, seqEnd(genome[v], c) - coord);
  decodeSpliceSignal(out, signal, alphabet.decode, alphabet.complement,
		     isSenseStrand == isForwardStrand);
}

void SplitAligner::spliceEndSignal(char *out,
				   unsigned alnNum, unsigned queryPos,
				   bool isSenseStrand) const {
  const UnsplitAlignment& a = alns[alnNum];
  bool isForwardStrand = a.isForwardStrand();
  StringNumMap::const_iterator f = chromosomeIndex.find(a.rname);
  size_t v = f->second % maxGenomeVolumes();
  size_t c = f->second / maxGenomeVolumes();
  uchar signal[2];
  unsigned coord = cell(spliceEndCoords, alnNum, queryPos);
  if (isForwardStrand) getPrevSignal(signal, seqBeg(genome[v], c) + coord);
  else                 getNextSignal(signal, seqEnd(genome[v], c) - coord);
  decodeSpliceSignal(out, signal, alphabet.decode, alphabet.complement,
		     isSenseStrand == isForwardStrand);
}

struct RnameAndStrandLess {
  RnameAndStrandLess(const UnsplitAlignment *a) : alns(a) {}

  bool operator()(unsigned a, unsigned b) const {
    return
      alns[a].qstrand != alns[b].qstrand ? alns[a].qstrand < alns[b].qstrand :
      std::strcmp(alns[a].rname, alns[b].rname) < 0;
  }

  const UnsplitAlignment *alns;
};

void SplitAligner::initRnameAndStrandIds() {
  rnameAndStrandIds.resize(numAlns);
  RnameAndStrandLess less(&alns[0]);
  stable_sort(sortedAlnIndices.begin(), sortedAlnIndices.end(), less);
  unsigned c = 0;
  for (unsigned i = 0; i < numAlns; ++i) {
    unsigned k = sortedAlnIndices[i];
    if (i > 0 && less(sortedAlnIndices[i-1], k)) ++c;
    rnameAndStrandIds[k] = c;
  }
}

void SplitAligner::dpExtensionMinScores(int maxJumpScore,
					size_t& minScore1,
					size_t& minScore2) const {
  if (!chromosomeIndex.empty()) maxJumpScore += maxSpliceBegEndScore;
  assert(maxJumpScore + insExistenceScore <= 0);
  minScore1 = 1 - (maxJumpScore + insExistenceScore);
  minScore2 = 1 - (maxJumpScore + maxJumpScore + insExistenceScore);
}

static size_t dpExtension(size_t maxScore, size_t minScore, size_t divisor) {
  return (maxScore > minScore) ? (maxScore - minScore) / divisor : 0;
}

void SplitAligner::initDpBounds() {
  minBeg = -1;
  for (unsigned i = 0; i < numAlns; ++i)
    minBeg = std::min(minBeg, alns[i].qstart);

  maxEnd = 0;
  for (unsigned i = 0; i < numAlns; ++i)
    maxEnd = std::max(maxEnd, alns[i].qend);

  dpBegs.resize(numAlns);
  dpEnds.resize(numAlns);

  // We will do dynamic programming along the length of each candidate
  // alignment.  But sometimes we need to consider "end gaps" and
  // extend the DP beyond the ends of each candidate.  Here we define
  // extensions, which aim to be as short as possible, but guarantee
  // to find the optimal split alignment score.  (Currently, they are
  // not as short as possible: this could be improved.)  We use these
  // facts:

  // The highest possible score for a given length is
  // length * maxMatchScore

  // An extension of length x must have a (negative) score <=
  // maxJumpScore + insExistenceScore + insExtensionScore * x

  assert(insExtensionScore < 0);
  assert(maxMatchScore >= 0);

  size_t oldDiv = -insExtensionScore;
  size_t newDiv = maxMatchScore - insExtensionScore;

  size_t minScore1 = -1;
  size_t minScore2 = -1;
  if (jumpProb > 0.0 || splicePrior > 0.0) {
    int m = (splicePrior > 0.0) ? maxSpliceScore : jumpScore;
    dpExtensionMinScores(m, minScore1, minScore2);
  }

  for (unsigned i = 0; i < numAlns; ++i) {
    size_t b = alns[i].qstart;
    size_t e = alns[i].qend;

    size_t bo = dpExtension(maxMatchScore * (e - b), minScore1, oldDiv);
    size_t bj = dpExtension(maxMatchScore * (maxEnd - b), minScore2, oldDiv);
    size_t bn = dpExtension(maxMatchScore * (b - minBeg), minScore1, newDiv);
    dpBegs[i] = b - std::min(std::max(bo, bj), bn);

    size_t eo = dpExtension(maxMatchScore * (e - b), minScore1, oldDiv);
    size_t ej = dpExtension(maxMatchScore * (e - minBeg), minScore2, oldDiv);
    size_t en = dpExtension(maxMatchScore * (maxEnd - e), minScore1, newDiv);
    dpEnds[i] = e + std::min(std::max(eo, ej), en);
  }

  // This sets the coordinate system for a ragged matrix, with numAlns
  // rows, where row i has cells from dpBeg(i) to dpEnd(i) inclusive.
  // (The final cell per row is used in some matrices but not others.)
  matrixRowOrigins.resize(numAlns);
  size_t s = 0;
  for (unsigned i = 0; i < numAlns; ++i) {
    s -= dpBeg(i);
    matrixRowOrigins[i] = s;
    s += dpEnd(i) + 1;
  }
}

void SplitAligner::layout(std::vector<UnsplitAlignment>::const_iterator beg,
			  std::vector<UnsplitAlignment>::const_iterator end) {
    assert(end > beg);
    numAlns = end - beg;
    alns = beg;

    sortedAlnIndices.resize(numAlns);
    for (unsigned i = 0; i < numAlns; ++i) sortedAlnIndices[i] = i;
    oldInplayAlnIndices.resize(numAlns);
    newInplayAlnIndices.resize(numAlns);

    rBegs.resize(numAlns);
    rEnds.resize(numAlns);

    if (splicePrior > 0.0 || !chromosomeIndex.empty()) {
	initRbegsAndEnds();
	//initRnameAndStrandIds();
    }
    initRnameAndStrandIds();

    initDpBounds();
}

size_t SplitAligner::memory(bool isViterbi, bool isBothSpliceStrands) const {
  size_t numOfStrands = isBothSpliceStrands ? 2 : 1;
  size_t x = 2 * sizeof(int) + 2 * sizeof(double);
  if (splicePrior > 0 || !chromosomeIndex.empty()) x += 2 * sizeof(unsigned);
  if (!chromosomeIndex.empty()) x += 2;
  if (isViterbi) x += sizeof(long) * numOfStrands;
  x += 2 * sizeof(double) * numOfStrands;
  return x * cellsPerDpMatrix();
}

void SplitAligner::initMatricesForOneQuery() {
  resizeMatrix(Amat);
  resizeMatrix(Dmat);
  resizeMatrix(Aexp);
  resizeMatrix(Dexp);

  for (unsigned i = 0; i < numAlns; i++) calcBaseScores(i);

  if (splicePrior > 0.0 || !chromosomeIndex.empty()) {
    resizeMatrix(spliceBegCoords);
    resizeMatrix(spliceEndCoords);
    for (unsigned i = 0; i < numAlns; ++i) initSpliceCoords(i);
  }

  if (!chromosomeIndex.empty()) {
    resizeMatrix(spliceBegSignals);
    resizeMatrix(spliceEndSignals);
    for (unsigned i = 0; i < numAlns; ++i) initSpliceSignals(i);
  }

  transform(Amat.begin(), Amat.end(), Aexp.begin(), scaledExp);
  transform(Dmat.begin(), Dmat.end(), Dexp.begin(), scaledExp);
  // if x/scale < about -745, then exp(x/scale) will be exactly 0.0
}

void SplitAligner::flipSpliceSignals() {
  Vmat.swap(VmatRev);
  Vvec.swap(VvecRev);
  Fmat.swap(FmatRev);
  Bmat.swap(BmatRev);
  rescales.swap(rescalesRev);

  for (int i = 0; i < 16; ++i) {
    int j = 15 - ((i%4) * 4 + (i/4));  // reverse-complement
    std::swap(spliceBegScores[i], spliceEndScores[j]);
    std::swap(spliceBegProbs[i], spliceEndProbs[j]);
  }
}

double SplitAligner::spliceSignalStrandLogOdds() const {
  assert(rescales.size() == rescalesRev.size());
  double logOdds = 0;
  for (unsigned j = 0; j < rescales.size(); ++j) {
    logOdds += std::log(rescales[j] / rescalesRev[j]);
  }
  return logOdds;
}

// 1st 1 million reads from SRR359290.fastq:
// lastal -Q1 -e120 hg19/last/female-1111110m
// last-split-probs -s150 -b.01 splicePrior=0
// distance sample size: 41829
// distance quartiles: 312 1122 3310
// estimated mean ln[distance] 7.02287
// estimated standard deviation of ln[distance] 1.75073
// This log-normal fits the data pretty well, especially for longer
// introns, but it's a bit inaccurate for short introns.

// last-split-probs -s150 splicePrior=0.01 meanLogDist=7.0 sdevLogDist=1.75
// distance sample size: 46107
// distance quartiles: 316 1108 3228
// estimated mean ln[distance] 7.01031
// estimated standard deviation of ln[distance] 1.72269

void SplitAligner::setSpliceParams(double splicePriorIn,
				   double meanLogDistIn,
				   double sdevLogDistIn) {
  splicePrior = splicePriorIn;
  meanLogDist = meanLogDistIn;
  sdevLogDist = sdevLogDistIn;

  if (splicePrior <= 0.0) return;

  const double rootTwoPi = std::sqrt(8.0 * std::atan(1.0));
  double s2 = sdevLogDist * sdevLogDist;
  spliceTerm1 = -std::log(sdevLogDist * rootTwoPi / splicePrior);
  spliceTerm2 = -0.5 / s2;

  double max1 = spliceTerm1 - meanLogDist + s2 * 0.5;
  int max2 = std::floor(scale * max1 + 0.5);
  maxSpliceScore = std::max(max2, jumpScore);

  // Set maxSpliceDist so as to ignore splices whose score would be
  // less than jumpScore.  By solving this quadratic equation:
  // spliceTerm1 + spliceTerm2 * (logDist - meanLogDist)^2 - logDist =
  // jumpScore / scale
  double r = s2 + 2 * (spliceTerm1 - meanLogDist - jumpScore / scale);
  if (r < 0) {
    maxSpliceDist = 0;
  } else {
    double logMode = meanLogDist - s2;  // ln(mode of log-normal distribution)
    double maxLogDist = logMode + sdevLogDist * std::sqrt(r);
    double maxDist = std::exp(maxLogDist);
    maxSpliceDist = -1;  // maximum possible unsigned value
    if (maxDist < maxSpliceDist) maxSpliceDist = std::floor(maxDist);
  }

  spliceTableSize = 256 * 256 * 64;
  spliceTableSize = std::min(spliceTableSize, maxSpliceDist);
  spliceScoreTable.resize(spliceTableSize);
  spliceProbTable.resize(spliceTableSize);
  for (unsigned i = 1; i < spliceTableSize; ++i) {
    int s = calcSpliceScore(i);
    spliceScoreTable[i] = s;
    spliceProbTable[i] = scaledExp(s);
  }
}

void SplitAligner::setParams(int gapExistenceScoreIn, int gapExtensionScoreIn,
			     int insExistenceScoreIn, int insExtensionScoreIn,
			     int jumpScoreIn, int restartScoreIn,
			     double scaleIn, int qualityOffsetIn) {
  gapExistenceScore = gapExistenceScoreIn;
  gapExtensionScore = gapExtensionScoreIn;
  insExistenceScore = insExistenceScoreIn;
  insExtensionScore = insExtensionScoreIn;
  jumpScore = jumpScoreIn;
  restartScore = restartScoreIn;
  scale = scaleIn;
  scaledExp.setBase(std::exp(1.0 / scale));
  qualityOffset = qualityOffsetIn;
  jumpProb = scaledExp(jumpScore);
  restartProb = scaledExp(restartScore);
}

static int scoreFromProb(double prob, double scale) {
  return std::floor(scale * std::log(prob) + 0.5);
}

void SplitAligner::setSpliceSignals() {
  // If an RNA-DNA alignment reaches position i in the DNA, the
  // probability of splicing from i to j is:
  //   P(i & j)  =  d(i) * a(j) * f(j - i),
  // where:
  // d(i) and a(j) depend on the DNA sequences at i and j, e.g. GT-AG,
  // and:
  // f(j - i) is a probability density function, e.g. log-normal.
  // So: the sum over j of f(j - i)  =  1.
  // The probability of splicing from i to anywhere is:
  //   P(i)  =  d(i) * sum over j of [a(j) * f(j - i)]
  // So, a typical value of P(i) is: typical(d) * typical(a).

  // Here, we set the values of d(i) and a(j).
  // XXX We should allow the user to choose different values.
  // Only the relative values matter, because we will normalize them
  // (so that the overall splice probability is set by splicePrior).

  // The values for non-GT-AG signals are unnaturally high, to allow
  // for various kinds of error.

  double dGT = 0.95;
  double dGC = 0.02;
  double dAT = 0.004;
  double dNN = 0.002;

  double aAG = 0.968;
  double aAC = 0.004;
  double aNN = 0.002;

  // We assume the dinucleotides have roughly equal 1/16 abundances.

  double dAvg = (dGT + dGC + dAT + dNN * 13) / 16;
  double aAvg = (aAG + aAC + aNN * 14) / 16;

  for (int i = 0; i < 17; ++i) {
    spliceBegScores[i] = scoreFromProb(dNN / dAvg, scale);
    spliceEndScores[i] = scoreFromProb(aNN / aAvg, scale);
  }

  spliceBegScores[2 * 4 + 3] = scoreFromProb(dGT / dAvg, scale);
  spliceBegScores[2 * 4 + 1] = scoreFromProb(dGC / dAvg, scale);
  spliceBegScores[0 * 4 + 3] = scoreFromProb(dAT / dAvg, scale);

  spliceEndScores[0 * 4 + 2] = scoreFromProb(aAG / aAvg, scale);
  spliceEndScores[0 * 4 + 1] = scoreFromProb(aAC / aAvg, scale);

  for (int i = 0; i < 17; ++i) {
    spliceBegProbs[i] = scaledExp(spliceBegScores[i]);
    spliceEndProbs[i] = scaledExp(spliceEndScores[i]);
  }

  maxSpliceBegEndScore = arrayMax(spliceBegScores) + arrayMax(spliceEndScores);
}

void SplitAligner::printParameters() const {
  if (jumpProb > 0.0) {
    std::cout << "# trans=" << jumpScore << "\n";
  }

  if (splicePrior > 0.0 && jumpProb > 0.0) {
    std::cout << "# cismax=" << maxSpliceDist << "\n";
  }

  if (!chromosomeIndex.empty()) {
    std::cout << "#"
	      << " GT=" << spliceBegScores[2 * 4 + 3]
	      << " GC=" << spliceBegScores[2 * 4 + 1]
	      << " AT=" << spliceBegScores[0 * 4 + 3]
	      << " NN=" << spliceBegScores[0 * 4 + 0]
	      << "\n";

    std::cout << "#"
	      << " AG=" << spliceEndScores[0 * 4 + 2]
	      << " AC=" << spliceEndScores[0 * 4 + 1]
	      << " NN=" << spliceEndScores[0 * 4 + 0]
	      << "\n";
  }
}

static void readPrjFile(const std::string& baseName,
			std::string& alphabetLetters,
			size_t& seqCount,
			size_t& volumes) {
  size_t fileBitsPerInt = 32;
  seqCount = volumes = -1;

  std::string fileName = baseName + ".prj";
  std::ifstream f(fileName.c_str());
  if (!f) err("can't open file: " + fileName);

  std::string line, word;
  while (getline(f, line)) {
    std::istringstream iss(line);
    getline(iss, word, '=');
    if (word == "alphabet") iss >> alphabetLetters;
    if (word == "numofsequences") iss >> seqCount;
    if (word == "volumes") iss >> volumes;
    if (word == "integersize") iss >> fileBitsPerInt;
  }

  if (alphabetLetters != "ACGT") err("can't read file: " + fileName);

  if (fileBitsPerInt != sizeof(MultiSequence::indexT) * CHAR_BIT) {
    if (fileBitsPerInt == 32) err("please use last-split for " + baseName);
    if (fileBitsPerInt == 64) err("please use last-split8 for " + baseName);
    err("weird integersize in " + fileName);
  }
}

void SplitAligner::readGenomeVolume(const std::string& baseName,
				    size_t seqCount,
				    size_t volumeNumber) {
  if (seqCount + 1 == 0) err("can't read: " + baseName);

  genome[volumeNumber].fromFiles(baseName, seqCount, 0);

  for (unsigned long long i = 0; i < seqCount; ++i) {
    char s = genome[volumeNumber].strand(i);
    if (s == '-') continue;
    std::string n = genome[volumeNumber].seqName(i);
    unsigned long long j = i * maxGenomeVolumes() + volumeNumber;
    if (!chromosomeIndex.insert(std::make_pair(n, j)).second)
      err("duplicate sequence name: " + n);
  }
}

void SplitAligner::readGenome(const std::string& baseName) {
  std::string alphabetLetters;
  size_t seqCount, volumes;
  readPrjFile(baseName, alphabetLetters, seqCount, volumes);

  if (volumes + 1 > 0 && volumes > 1) {
    if (volumes > maxGenomeVolumes()) err("too many volumes: " + baseName);
    for (size_t i = 0; i < volumes; ++i) {
      std::string b = baseName + stringify(i);
      size_t c, v;
      readPrjFile(b, alphabetLetters, c, v);
      readGenomeVolume(b, c, i);
    }
  } else {
    readGenomeVolume(baseName, seqCount, 0);
  }

  alphabet.fromString(alphabetLetters);
}

static double probFromPhred(double s) {
  return std::pow(10.0, -0.1 * s);
}

static int generalizedScore(double score, double scale, double phredScore,
			    double letterProb) {
  double r = std::exp(score / scale);
  double p = probFromPhred(phredScore);
  if (p >= 1) p = 0.999999;  // kludge to avoid numerical instability
  double otherProb = 1 - letterProb;
  assert(otherProb > 0);
  double u = p / otherProb;
  double x = (1 - u) * r + u;
  assert(x > 0);
  return std::floor(scale * std::log(x) + 0.5);
}

static int max(const std::vector< std::vector<int> >& matrix) {
  int m = matrix.at(0).at(0);
  for (unsigned i = 0; i < matrix.size(); ++i)
    for (unsigned j = 0; j < matrix[i].size(); ++j)
      m = std::max(m, matrix[i][j]);
  return m;
}

static int min(const std::vector< std::vector<int> >& matrix) {
  int m = matrix.at(0).at(0);
  for (unsigned i = 0; i < matrix.size(); ++i)
    for (unsigned j = 0; j < matrix[i].size(); ++j)
      m = std::min(m, matrix[i][j]);
  return m;
}

static int matrixLookup(const std::vector< std::vector<int> >& matrix,
			const std::string& rowNames,
			const std::string& colNames, char x, char y) {
  std::string::size_type row = rowNames.find(x);
  std::string::size_type col = colNames.find(y);
  if (row == std::string::npos || col == std::string::npos)
    return min(matrix);
  else
    return matrix.at(row).at(col);
}

void SplitAligner::setScoreMat(const std::vector< std::vector<int> >& matrix,
			       const std::string& rowNames,
			       const std::string& colNames) {
  const std::string bases = "ACGT";

  // Reverse-engineer the abundances of ACGT from the score matrix:
  unsigned blen = bases.size();
  std::vector<int> bvec(blen * blen);
  std::vector<int *> bmat(blen);
  for (unsigned i = 0; i < blen; ++i) bmat[i] = &bvec[i * blen];
  for (unsigned i = 0; i < blen; ++i)
    for (unsigned j = 0; j < blen; ++j)
      bmat[i][j] = matrixLookup(matrix, rowNames, colNames,
				bases[i], bases[j]);

  mcf::SubstitutionMatrixStats stats;
  stats.calcFromScale(&bmat[0], blen, scale);

  for (int i = 64; i < 128; ++i) {
    char x = std::toupper(i);
    for (int j = 64; j < 128; ++j) {
      char y = std::toupper(j);
      int score = matrixLookup(matrix, rowNames, colNames, x, y);
      for (int q = 0; q < numQualCodes; ++q) {
	std::string::size_type xc = bases.find(x);
	std::string::size_type yc = bases.find(y);
	if (xc == std::string::npos || yc == std::string::npos) {
	  score_mat[i % 64][j % 64][q] = score;
	} else {
	  double p = stats.letterProbs2()[yc];
	  score_mat[i % 64][j % 64][q] = generalizedScore(score, scale, q, p);
	}
      }
    }
  }

  maxMatchScore = max(matrix);
}

}
