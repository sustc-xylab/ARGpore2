// Copyright 2010, 2013 Martin C. Frith

#include "gaplessPssmXdrop.hh"
#include <stdexcept>

static void err(const char *s) { throw std::overflow_error(s); }

namespace cbrc {

int forwardGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*pssm++)[*seq++];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in forward gapless extension with PSSM");
  return score;
}

int reverseGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*--pssm)[*--seq];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in reverse gapless extension with PSSM");
  return score;
}

const uchar *forwardGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score) {
  int s = 0;
  while (s < score) s += (*pssm++)[*seq++];
  return seq;
}

const uchar *reverseGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score) {
  int s = 0;
  while (s < score) s += (*--pssm)[*--seq];
  return seq;
}

bool isOptimalGaplessPssmXdrop(const uchar *seq,
                               const uchar *seqEnd,
                               const ScoreMatrixRow *pssm,
                               int maxScoreDrop) {
  int score = 0;
  int maxScore = 0;
  while (seq < seqEnd) {
    score += (*pssm++)[*seq++];
    if (score > maxScore) maxScore = score;
    else if (score <= 0 ||                       // non-optimal prefix
             seq == seqEnd ||                    // non-optimal suffix
             score < maxScore - maxScoreDrop) {  // excessive score drop
      return false;
    }
  }
  return true;
}

int gaplessPssmXdropOverlap(const uchar *seq,
			    const ScoreMatrixRow *pssm,
			    int maxScoreDrop,
			    size_t &reverseLength,
			    size_t &forwardLength) {
  int minScore = 0;
  int maxScore = 0;
  int score = 0;

  const uchar *rs = seq;
  const ScoreMatrixRow *rp = pssm;
  while (true) {
    --rs;  --rp;
    int s = (*rp)[*rs];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    else if (score < minScore) minScore = score;
  }

  maxScore = score - minScore;

  const uchar *fs = seq;
  const ScoreMatrixRow *fp = pssm;
  while (true) {
    int s = (*fp)[*fs];
    if (s <= -INF) break;
    score += s;
    if (score > maxScore) maxScore = score;
    else if (score < maxScore - maxScoreDrop) return -INF;
    ++fs;  ++fp;
  }

  reverseLength = seq - (rs + 1);
  forwardLength = fs - seq;
  return score;
}

int gaplessPssmAlignmentScore(const uchar *seq,
                              const uchar *seqEnd,
                              const ScoreMatrixRow *pssm) {
  int score = 0;
  while (seq < seqEnd) score += (*pssm++)[*seq++];
  return score;
}

}
