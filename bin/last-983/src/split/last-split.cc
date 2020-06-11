// Copyright 2013, 2014 Martin C. Frith

#include "last-split.hh"

#include "cbrc_split_aligner.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstring>  // strcmp
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static std::istream& openIn(const std::string& fileName, std::ifstream& ifs) {
  if (fileName == "-") return std::cin;
  ifs.open(fileName.c_str());
  if (!ifs) err("can't open file: " + fileName);
  return ifs;
}

// Does the string start with the prefix?
static bool startsWith(const std::string& s, const char* prefix) {
  const char* t = s.c_str();
  for (;;) {
    if (*prefix == 0) return true;
    if (*prefix != *t) return false;
    ++t;
    ++prefix;
  }
}

// Does the string have no non-space characters?
static bool isSpace(const std::string& s) {
  const char* t = s.c_str();
  for (;;) {
    if (*t == 0) return true;
    if (!std::isspace(*t)) return false;
    ++t;
  }
}

static int scoreFromProb(double prob, double scale) {
  return std::floor(scale * std::log(prob) + 0.5);
}

// Defines an ordering, for sorting.
static bool less(const cbrc::UnsplitAlignment& a,
		 const cbrc::UnsplitAlignment& b) {
  int qnameCmp = std::strcmp(a.qname, b.qname);
  if (qnameCmp  != 0        ) return qnameCmp  < 0;
  if (a.qstart  != b.qstart ) return a.qstart  < b.qstart;
  if (a.qend    != b.qend   ) return a.qend    < b.qend;
  if (a.qstrand != b.qstrand) return a.qstrand < b.qstrand;
  int qalignCmp = std::strcmp(a.qalign, b.qalign);
  if (qalignCmp != 0        ) return qalignCmp < 0;
  int ralignCmp = std::strcmp(a.ralign, b.ralign);
  if (ralignCmp != 0        ) return ralignCmp < 0;
  return a.linesBeg < b.linesBeg;  // stabilizes the sort
}

static void printSense(double senseStrandLogOdds) {
  double b = senseStrandLogOdds / std::log(2.0);
  if (b < 0.1 && b > -0.1) b = 0;
  else if (b > 10) b = std::floor(b + 0.5);
  else if (b < -10) b = std::ceil(b - 0.5);
  int precision = (b < 10 && b > -10) ? 2 : 3;
  std::cout << std::setprecision(precision) << " sense=" << b;
}

static void doOneAlignmentPart(cbrc::SplitAligner& sa,
			       const cbrc::UnsplitAlignment& a,
			       unsigned numOfParts, unsigned partNum,
			       unsigned alnNum,
			       unsigned qSliceBeg, unsigned qSliceEnd,
			       bool isSenseStrand, double senseStrandLogOdds,
			       const LastSplitOptions& opts,
			       bool isAlreadySplit) {
  if (qSliceBeg >= a.qend || qSliceEnd <= a.qstart) {
    return;  // this can happen for spliced alignment!
  }

  unsigned qSliceBegTrimmed = qSliceBeg;
  unsigned qSliceEndTrimmed = qSliceEnd;
  unsigned alnBeg, alnEnd;
  cbrc::mafSliceBeg(a.ralign, a.qalign, a.qstart, qSliceBegTrimmed, alnBeg);
  cbrc::mafSliceEnd(a.ralign, a.qalign, a.qend,   qSliceEndTrimmed, alnEnd);

  if (qSliceBegTrimmed >= qSliceEndTrimmed) {
    return;  // I think this can happen for spliced alignment
  }

  int score =
    sa.segmentScore(alnNum, qSliceBeg, qSliceEnd) -
    sa.segmentScore(alnNum, qSliceBeg, qSliceBegTrimmed) -
    sa.segmentScore(alnNum, qSliceEndTrimmed, qSliceEnd);
  if (score < opts.score) return;

  std::vector<double> p;
  if (opts.direction != 0) {
    p = sa.marginalProbs(qSliceBegTrimmed, alnNum, alnBeg, alnEnd);
  }
  std::vector<double> pRev;
  if (opts.direction != 1) {
    sa.flipSpliceSignals();
    pRev = sa.marginalProbs(qSliceBegTrimmed, alnNum, alnBeg, alnEnd);
    sa.flipSpliceSignals();
  }
  if (opts.direction == 0) p.swap(pRev);
  if (opts.direction == 2) {
    double reverseProb = 1 / (1 + std::exp(senseStrandLogOdds));
    // the exp might overflow to inf, but that should be OK
    double forwardProb = 1 - reverseProb;
    for (unsigned i = 0; i < p.size(); ++i) {
      p[i] = forwardProb * p[i] + reverseProb * pRev[i];
    }
  }

  assert(!p.empty());
  double mismap = 1.0 - *std::max_element(p.begin(), p.end());
  mismap = std::max(mismap, 1e-10);
  if (mismap > opts.mismap) return;
  int mismapPrecision = 3;

  std::vector<std::string> s = cbrc::mafSlice(a.linesBeg, a.linesEnd,
					      alnBeg, alnEnd);
  s.push_back(cbrc::pLineFromProbs(p));

  if (isAlreadySplit && s.end()[-2][0] == 'p') {
    mismap = cbrc::pLinesToErrorProb(s.end()[-2].c_str(), s.end()[-1].c_str());
    if (mismap > opts.mismap) return;
    mismapPrecision = 2;
    if (opts.format == 'm') s.pop_back();
  }

  if (opts.format == 'm') s.pop_back();

  if (opts.no_split && a.linesBeg[0][0] == 'a') {
    std::cout << a.linesBeg[0];
  } else {
    std::cout << "a score=" << score;
  }
  std::cout << std::setprecision(mismapPrecision) << " mismap=" << mismap;
  if (opts.direction == 2) printSense(senseStrandLogOdds);
  if (!opts.genome.empty() && !opts.no_split) {
    char signal[3] = {0};
    if (partNum > 0) {
      sa.spliceEndSignal(signal, alnNum, qSliceBeg, isSenseStrand);
      std::cout << (isSenseStrand ? " acc=" : " don=") << signal;
    }
    if (partNum + 1 < numOfParts) {
      sa.spliceBegSignal(signal, alnNum, qSliceEnd, isSenseStrand);
      std::cout << (isSenseStrand ? " don=" : " acc=") << signal;
    }
  }
  std::cout << "\n" << std::setprecision(6);

  if (a.isFlipped()) cbrc::flipMafStrands(s.begin(), s.end());
  if (opts.no_split && a.linesEnd[-1][0] == 'c') s.push_back(a.linesEnd[-1]);
  cbrc::printMaf(s);
}

static void doOneQuery(std::vector<cbrc::UnsplitAlignment>::const_iterator beg,
		       std::vector<cbrc::UnsplitAlignment>::const_iterator end,
		       cbrc::SplitAligner& sa, const LastSplitOptions& opts,
		       bool isAlreadySplit) {
  if (opts.verbose) std::cerr << beg->qname << "\t" << (end - beg);
  sa.layout(beg, end);
  if (opts.verbose) std::cerr << "\tcells=" << sa.cellsPerDpMatrix();
  size_t bytes = sa.memory(!opts.no_split, opts.direction == 2);
  if (bytes > opts.bytes) {
    if (opts.verbose) std::cerr << "\n";
    std::cerr << "last-split: skipping sequence " << beg->qname
	      << " (" << bytes << " bytes)\n";
    return;
  }
  sa.initMatricesForOneQuery();

  if (opts.direction != 0) {
    sa.forward();
    sa.backward();
  }
  if (opts.direction != 1) {
    sa.flipSpliceSignals();
    sa.forward();
    sa.backward();
    sa.flipSpliceSignals();
  }

  double senseStrandLogOdds =
    (opts.direction == 2) ? sa.spliceSignalStrandLogOdds() : 0;

  if (opts.no_split) {
    if (opts.verbose) std::cerr << "\n";
    unsigned numOfParts = end - beg;
    for (unsigned i = 0; i < numOfParts; ++i) {
      doOneAlignmentPart(sa, beg[i], numOfParts, i, i,
			 beg[i].qstart, beg[i].qend,
			 true, senseStrandLogOdds, opts, isAlreadySplit);
    }
  } else {
    long viterbiScore = LONG_MIN;
    if (opts.direction != 0) {
      viterbiScore = sa.viterbi();
      if (opts.verbose) std::cerr << "\t" << viterbiScore;
    }
    long viterbiScoreRev = LONG_MIN;
    if (opts.direction != 1) {
      sa.flipSpliceSignals();
      viterbiScoreRev = sa.viterbi();
      sa.flipSpliceSignals();
      if (opts.verbose) std::cerr << "\t" << viterbiScoreRev;
    }
    bool isSenseStrand = (viterbiScore >= viterbiScoreRev);
    std::vector<unsigned> alnNums;
    std::vector<unsigned> queryBegs;
    std::vector<unsigned> queryEnds;
    if (isSenseStrand) {
      sa.traceBack(viterbiScore, alnNums, queryBegs, queryEnds);
    } else {
      sa.flipSpliceSignals();
      sa.traceBack(viterbiScoreRev, alnNums, queryBegs, queryEnds);
      sa.flipSpliceSignals();
    }
    std::reverse(alnNums.begin(), alnNums.end());
    std::reverse(queryBegs.begin(), queryBegs.end());
    std::reverse(queryEnds.begin(), queryEnds.end());

    if (opts.verbose) std::cerr << "\n";
    unsigned numOfParts = alnNums.size();
    for (unsigned k = 0; k < numOfParts; ++k) {
      unsigned i = alnNums[k];
      doOneAlignmentPart(sa, beg[i], numOfParts, k, i,
			 queryBegs[k], queryEnds[k],
			 isSenseStrand, senseStrandLogOdds,
			 opts, isAlreadySplit);
    }
  }
}

static void doOneBatch(std::vector<std::string>& mafLines,
		       const std::vector<unsigned>& mafEnds,
                       cbrc::SplitAligner& sa, const LastSplitOptions& opts,
		       bool isAlreadySplit) {
  std::vector<cbrc::UnsplitAlignment> mafs;
  mafs.reserve(mafEnds.size() - 1);  // saves memory: no excess capacity
  for (unsigned i = 1; i < mafEnds.size(); ++i)
    mafs.push_back(cbrc::UnsplitAlignment(mafLines.begin() + mafEnds[i-1],
					  mafLines.begin() + mafEnds[i]));

  sort(mafs.begin(), mafs.end(), less);
  std::vector<cbrc::UnsplitAlignment>::const_iterator b = mafs.begin();
  std::vector<cbrc::UnsplitAlignment>::const_iterator e = mafs.begin();
  size_t qendMax = 0;
  while (e < mafs.end()) {
    if (e->qend > qendMax) qendMax = e->qend;
    ++e;
    if (e == mafs.end() || std::strcmp(e->qname, b->qname) != 0 ||
	(e->qstart >= qendMax && !opts.isSplicedAlignment)) {
      doOneQuery(b, e, sa, opts, isAlreadySplit);
      b = e;
      qendMax = 0;
    }
  }
}

static void printParameters(const LastSplitOptions& opts) {
  std::cout << std::setprecision(12) << "#"
	    << " m=" << opts.mismap
	    << " s=" << opts.score;
  if (opts.isSplicedAlignment) {
    std::cout << " d=" << opts.direction
	      << " c=" << opts.cis
	      << " t=" << opts.trans
	      << " M=" << opts.mean
	      << " S=" << opts.sdev;
  }
  std::cout << "\n" << std::setprecision(6);
}

static void addMaf(std::vector<unsigned>& mafEnds,
		   const std::vector<std::string>& mafLines) {
  if (mafLines.size() > mafEnds.back())  // if we have new maf lines:
    mafEnds.push_back(mafLines.size());  // store the new end
}

void lastSplit(LastSplitOptions& opts) {
  cbrc::SplitAligner sa;
  std::vector< std::vector<int> > scoreMatrix;
  std::string rowNames, colNames;
  std::string line, word, name, key;
  int state = 0;
  int sequenceFormat = -1;
  int gapExistenceCost = -1;
  int gapExtensionCost = -1;
  int insExistenceCost = -1;
  int insExtensionCost = -1;
  int lastalScoreThreshold = -1;
  double scale = 0;
  double genomeSize = 0;
  std::vector<std::string> mafLines;  // lines of multiple MAF blocks
  std::vector<unsigned> mafEnds(1);  // which lines are in which MAF block
  bool isAlreadySplit = false;  // has the input already undergone last-split?

  for (unsigned i = 0; i < opts.inputFileNames.size(); ++i) {
    std::ifstream inFileStream;
    std::istream& input = openIn(opts.inputFileNames[i], inFileStream);
    while (getline(input, line)) {
      if (state == -1) {  // we are reading the score matrix within the header
	std::istringstream ls(line);
	std::vector<int> row;
	int score;
	ls >> word >> name;
	while (ls >> score) row.push_back(score);
	if (word == "#" && name.size() == 1 && !row.empty() && ls.eof()) {
	  rowNames.push_back(std::toupper(name[0]));
	  scoreMatrix.push_back(row);
	} else {
	  state = 0;
	}
      }
      if (state == 0) {  // we are reading the header
	std::istringstream ls(line);
	std::string names;
	ls >> word;
	while (ls >> name) {
	  if (name.size() == 1) names.push_back(std::toupper(name[0]));
	  else break;
	}
	if (word == "#" && !names.empty() && !ls && scoreMatrix.empty()) {
	  colNames = names;
	  state = -1;
	} else if (startsWith(line, "#")) {
	  std::istringstream ls(line);
	  while (ls >> word) {
	    std::istringstream ws(word);
	    getline(ws, key, '=');
	    if (key == "a") ws >> gapExistenceCost;
	    if (key == "b") ws >> gapExtensionCost;
	    if (key == "A") ws >> insExistenceCost;
	    if (key == "B") ws >> insExtensionCost;
	    if (key == "e") ws >> lastalScoreThreshold;
	    if (key == "t") ws >> scale;
	    if (key == "Q") ws >> sequenceFormat;
	    if (key == "letters") ws >> genomeSize;
	  }
	  // try to determine if last-split was already run (fragile):
	  if (startsWith(line, "# m=")) isAlreadySplit = true;
	} else if (!isSpace(line)) {
	  if (scoreMatrix.empty())
	    err("I need a header with score parameters");
	  if (gapExistenceCost < 0 || gapExtensionCost < 0 ||
	      insExistenceCost < 0 || insExtensionCost < 0 ||
	      lastalScoreThreshold < 0 || scale <= 0 || genomeSize <= 0)
	    err("can't read the header");
	  if (sequenceFormat == 2 || sequenceFormat >= 4)
	    err("unsupported Q format");
	  if (opts.score < 0)
	    opts.score = lastalScoreThreshold +
	      (opts.isSplicedAlignment ? scoreFromProb(100, scale) : 0);
	  int restartCost =
	    opts.isSplicedAlignment ? -(INT_MIN/2) : opts.score - 1;
	  double jumpProb = opts.isSplicedAlignment
	    ? opts.trans / (2 * genomeSize)  // 2 strands
	    : 0.0;
	  int jumpCost =
	    (jumpProb > 0.0) ? -scoreFromProb(jumpProb, scale) : -(INT_MIN/2);
	  int qualityOffset = (sequenceFormat == 3) ? 64 : 33;
	  printParameters(opts);
	  sa.setParams(-gapExistenceCost, -gapExtensionCost,
		       -insExistenceCost, -insExtensionCost,
		       -jumpCost, -restartCost, scale, qualityOffset);
	  double splicePrior = opts.isSplicedAlignment ? opts.cis : 0.0;
	  sa.setSpliceParams(splicePrior, opts.mean, opts.sdev);
	  sa.setScoreMat(scoreMatrix, rowNames, colNames);
	  sa.setSpliceSignals();
	  if (!opts.genome.empty()) sa.readGenome(opts.genome);
	  sa.printParameters();
	  std::cout << "#\n";
	  state = 1;
	}
      }
      if (state == 1) {  // we are reading alignments
	if (startsWith(line, "# batch")) {
	  addMaf(mafEnds, mafLines);
	  doOneBatch(mafLines, mafEnds, sa, opts, isAlreadySplit);
	  mafLines.clear();
	  mafEnds.resize(1);
	} else if (isSpace(line)) {
	  addMaf(mafEnds, mafLines);
	} else if (std::strchr(opts.no_split ? "asqpc" : "sqpc", line[0])) {
	  mafLines.push_back(line);
	}
      }
      if (startsWith(line, "#") && !startsWith(line, "# batch"))
	std::cout << line << "\n";
    }
  }
  addMaf(mafEnds, mafLines);
  doOneBatch(mafLines, mafEnds, sa, opts, isAlreadySplit);
}
