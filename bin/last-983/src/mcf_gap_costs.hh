// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

// Gap costs for pair-wise sequence alignment.
// A gap of length k costs: openCost + k * growCost.

// "Generalized affine gap cost" is supported [SF Altschul 1998
// Proteins 32(1):88-96].  In this scheme, a "gap" may consist of
// unaligned regions in both sequences.  If these unaligned regions
// have sizes j and k, where j <= k, the cost is:
// openCost + growCost*(k-j) + pairCost*j

// Different costs for deletions and insertions are supported.

// Piecewise linear gap costs are supported.

#ifndef MCF_GAP_COSTS_HH
#define MCF_GAP_COSTS_HH

#include <vector>

namespace mcf {

struct GapCosts {
  struct Piece {
    int openCost;
    int growCost;
  };

  std::vector<Piece> delPieces;
  std::vector<Piece> insPieces;
  int pairCost;

  // Assign piecewise linear open and grow costs, and one pairCost.
  // If unalignedPairCost <= 0, assign non-generalized costs.
  // Throw a runtime_error if any growCost is <= 0.
  // delOpenCosts.size() must equal delGrowCosts.size(), and
  // insOpenCosts.size() must equal insGrowCosts.size().
  // The vectors must not be empty.
  void assign(const std::vector<int> &delOpenCosts,
	      const std::vector<int> &delGrowCosts,
	      const std::vector<int> &insOpenCosts,
	      const std::vector<int> &insGrowCosts,
	      int unalignedPairCost);

  // The cost of a "gap" consisting of unaligned letters in the query
  // and/or reference sequence
  int cost(int refInsertLen, int qryInsertLen) const;
};

}

#endif
