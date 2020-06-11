// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mcf_gap_costs.hh"
#include <stdexcept>
#include <assert.h>
#include <limits.h>
#include <stddef.h>  // size_t

static void err(const char *s) {
  throw std::runtime_error(s);
}

namespace mcf {

static void assignGapCostPieces(const std::vector<int> &openCosts,
				const std::vector<int> &growCosts,
				std::vector<GapCosts::Piece> &pieces) {
  assert(openCosts.size() == growCosts.size());
  pieces.clear();
  for (size_t i = 0; i < openCosts.size(); ++i) {
    GapCosts::Piece p = {openCosts[i], growCosts[i]};
    if (p.growCost <= 0) err("gap extension cost must be > 0");
    pieces.push_back(p);
  }
  assert(!pieces.empty());
}

void GapCosts::assign(const std::vector<int> &delOpenCosts,
		      const std::vector<int> &delGrowCosts,
		      const std::vector<int> &insOpenCosts,
		      const std::vector<int> &insGrowCosts,
		      int unalignedPairCost) {
  assignGapCostPieces(delOpenCosts, delGrowCosts, delPieces);
  assignGapCostPieces(insOpenCosts, insGrowCosts, insPieces);
  if (unalignedPairCost > 0) {
    pairCost = unalignedPairCost;
  } else {
    pairCost = delPieces[0].openCost + delPieces[0].growCost +
      insPieces[0].openCost + insPieces[0].growCost;
  }
}

int GapCosts::cost(int refInsertLen, int qryInsertLen) const {
  int genCost = INT_MAX;

  int delCost = 0;
  if (refInsertLen > 0) {
    delCost = INT_MAX;
    for (size_t i = 0; i < delPieces.size(); ++i) {
      int s = delPieces[i].openCost + delPieces[i].growCost * refInsertLen;
      delCost = std::min(delCost, s);
      if (refInsertLen >= qryInsertLen) {
	int t = s + (pairCost - delPieces[i].growCost) * qryInsertLen;
	genCost = std::min(genCost, t);
      }
    }
  }

  int insCost = 0;
  if (qryInsertLen > 0) {
    insCost = INT_MAX;
    for (size_t i = 0; i < insPieces.size(); ++i) {
      int s = insPieces[i].openCost + insPieces[i].growCost * qryInsertLen;
      insCost = std::min(insCost, s);
      if (qryInsertLen >= refInsertLen) {
	int t = s + (pairCost - insPieces[i].growCost) * refInsertLen;
	genCost = std::min(genCost, t);
      }
    }
  }

  return std::min(delCost + insCost, genCost);
}

}
