// Copyright 2016 Martin C. Frith

#include "SubsetMinimizerFinder.hh"
#include "CyclicSubsetSeed.hh"

namespace cbrc {

void SubsetMinimizerFinder::init(const CyclicSubsetSeed &seed,
				 const uchar *text,
				 size_t beg,
				 size_t end) {
  const uchar *m = seed.firstMap();
  while (beg < end && m[text[beg]] == CyclicSubsetSeed::DELIMITER) ++beg;
  minima.assign(1, beg);
}

bool SubsetMinimizerFinder::isMinimizer(const CyclicSubsetSeed &seed,
					const uchar *text,
					size_t pos,
					size_t end,
					size_t window) {
  const uchar *subsetMap = seed.firstMap();

  while (true) {
    size_t currentMinimum = minima[0];
    if (currentMinimum > pos) return false;
    if (currentMinimum == pos) return true;
    size_t newPos = minima.back() + 1;
    if (newPos == end) return false;
    const uchar *newText = text + newPos;
    if (subsetMap[*newText] == CyclicSubsetSeed::DELIMITER) {
      init(seed, text, newPos + 1, end);
      continue;
    }
    size_t stop = (newPos - currentMinimum >= window);
    size_t i = minima.size();
    while (i > stop) {
      size_t oldPos = minima[i - 1];
      const uchar *oldText = text + oldPos;
      if (!seed.isLess(newText, oldText, subsetMap)) break;
      --i;
    }
    minima.resize(i);
    if (stop) minima.erase(minima.begin());
    minima.push_back(newPos);
  }
}

}
