// Copyright 2008, 2009, 2010, 2011, 2013 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SubsetSuffixArray.hh"
#include <algorithm>  // iter_swap, min
//#include <iostream>  // for debugging

#ifdef HAS_CXX_THREADS
#include <thread>
#endif

using namespace cbrc;

namespace{
  typedef SubsetSuffixArray::indexT indexT;
  typedef SubsetSuffixArray::Range Range;
}

static void pushRange(std::vector<Range> &v,
		      indexT *beg, indexT *end, indexT depth) {
  if (end - beg > 1) {
    Range r = {beg, end, depth};
    v.push_back(r);
  }
}

static void insertionSort( const uchar* text, const CyclicSubsetSeed& seed,
			   indexT* beg, indexT* end, const uchar* subsetMap ){
  for( indexT* i = beg+1; i < end; ++i ){
    const uchar* newText = text + *i;
    for( indexT* j = i; j > beg; --j ){
      indexT* k = j - 1;
      const uchar* oldText = text + *k;
      if( !seed.isLess( newText, oldText, subsetMap ) ) break;
      std::iter_swap( j, k );
    }
  }
}

void SubsetSuffixArray::sort2( const uchar* text, indexT* beg,
			       const uchar* subsetMap ){
  indexT* mid = beg + 1;

  const uchar* s = text + *beg;
  const uchar* t = text + *mid;
  while( true ){
    uchar x = subsetMap[ *s ];
    uchar y = subsetMap[ *t ];
    if( x != y ){
      if( x > y ) std::iter_swap( beg, mid );
      break;
    }
    if( x == CyclicSubsetSeed::DELIMITER ) return;
    ++s;
    ++t;
    subsetMap = seed.nextMap(subsetMap);
  }

  if( isChildDirectionForward( beg ) ){
    setChildForward( beg + 0, mid );
  }else{
    setChildReverse( beg + 2, mid );
  }
}

// Specialized sort for 1 symbol + 1 delimiter.
// E.g. wildcard positions in spaced seeds.
void SubsetSuffixArray::radixSort1( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* begN = end;  // beginning of delimiters

  while( end0 < begN ){
    const indexT x = *end0;
    switch( subsetMap[ text[x] ] ){
    case 0:
      end0++;
      break;
    default:  // the delimiter subset
      *end0 = *--begN;
      *begN = x;
      break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
  }
}

// Specialized sort for 2 symbols + 1 delimiter.
// E.g. transition-constrained positions in subset seeds.
void SubsetSuffixArray::radixSort2( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* begN = end;  // beginning of delimiters

  while( end1 < begN ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      default:  // the delimiter subset
        *end1 = *--begN;
        *begN = x;
        break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 3 symbols + 1 delimiter.
// E.g. subset seeds for bisulfite-converted DNA.
void SubsetSuffixArray::radixSort3( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* beg2 = end;  // beginning of '2's
  indexT* begN = end;  // beginning of delimiters

  while( end1 < beg2 ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      case 2:
        *end1 = *--beg2;
        *beg2 = x;
        break;
      default:  // the delimiter subset
        *end1 = *--beg2;
        *beg2 = *--begN;
        *begN = x;
        break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, beg2, begN, depth );  // the '2's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( begN == end ) return;
    setChildForward( beg2, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg2 == beg ) return;
    setChildReverse( begN, beg2 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 4 symbols + 1 delimiter.  E.g. DNA.
void SubsetSuffixArray::radixSort4( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* end2 = beg;  // end of '2's
  indexT* beg3 = end;  // beginning of '3's
  indexT* begN = end;  // beginning of delimiters

  while( end2 < beg3 ){
    const indexT x = *end2;
    switch( subsetMap[ text[x] ] ){
    case 0:
      *end2++ = *end1;
      *end1++ = *end0;
      *end0++ = x;
      break;
    case 1:
      *end2++ = *end1;
      *end1++ = x;
      break;
    case 2:
      end2++;
      break;
    case 3:
      *end2 = *--beg3;
      *beg3 = x;
      break;
    default:  // the delimiter subset
      *end2 = *--beg3;
      *beg3 = *--begN;
      *begN = x;
      break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, end1, end2, depth );  // the '2's
  pushRange( rangeStack, beg3, begN, depth );  // the '3's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( end2 == end ) return;
    setChildForward( end1, end2 );
    if( begN == end ) return;
    setChildForward( beg3, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg3 == beg ) return;
    setChildReverse( begN, beg3 );
    if( end1 == beg ) return;
    setChildReverse( end2, end1 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

const unsigned numOfBuckets = 256;

void SubsetSuffixArray::radixSortN( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth,
				    unsigned subsetCount, indexT* bucketSize ){
  indexT* bucketEnd[numOfBuckets];

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for( indexT* i = beg; i < end; /* noop */ ){
    uchar oracle[256];
    uchar* oracleEnd =
      oracle + std::min( sizeof(oracle), size_t(end - i) );
    for( uchar* j = oracle; j < oracleEnd; ++j )
      *j = subsetMap[ text[ *i++ ] ];
    for( uchar* j = oracle; j < oracleEnd; ++j )
      ++bucketSize[ *j ];
  }

  // get bucket ends, and put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  indexT* pos = beg;
  for( unsigned i = 0; i < subsetCount; ++i ){
    indexT* nextPos = pos + bucketSize[i];
    pushRange( rangeStack, pos, nextPos, depth );
    pos = nextPos;
    bucketEnd[i] = pos;
  }
  // don't sort within the delimiter bucket:
  bucketEnd[ CyclicSubsetSeed::DELIMITER ] = end;

  if( isChildDirectionForward( beg ) ){
    pos = beg;
    for( unsigned i = 0; i < subsetCount; ++i ){
      indexT* nextPos = bucketEnd[i];
      if( nextPos == end ) break;
      setChildForward( pos, nextPos );
      pos = nextPos;
    }
  }else{
    pos = end;
    for( unsigned i = subsetCount; i > 0; --i ){
      indexT* nextPos = bucketEnd[i - 1];
      if( nextPos == beg ) break;
      setChildReverse( pos, nextPos );
      pos = nextPos;
    }
  }

  // permute items into the correct buckets:
  for( indexT* i = beg; i < end; /* noop */ ) {
    unsigned subset;  // unsigned is faster than uchar!
    indexT holdOut = *i;
    while( --bucketEnd[ subset = subsetMap[ text[holdOut] ] ] > i ){
      std::swap( *bucketEnd[subset], holdOut );
    }
    *i = holdOut;
    i += bucketSize[subset];
    bucketSize[subset] = 0;  // reset it so we can reuse it
  }
}

static size_t rangeSize(const Range &r) {
  return r.end - r.beg;
}

static size_t nextRangeSize(const std::vector<Range> &ranges) {
  return rangeSize(ranges.back());
}

static size_t rangeSizeSum(const std::vector<Range> &ranges) {
  size_t s = 0;
  for (size_t i = 0; i < ranges.size(); ++i) {
    s += rangeSize(ranges[i]);
  }
  return s;
}

static size_t numOfThreadsForOneRange(size_t numOfThreads,
				      size_t sizeOfThisRange,
				      size_t sizeOfAllRanges) {
  // We want:
  // min(t | sizeOfThisRange / t < sizeOfOtherRanges / (numOfThreads - (t+1)))
  // Or equivalently:
  // max(t | sizeOfThisRange / (t-1) >= sizeOfOtherRanges / (numOfThreads - t))
  double x = numOfThreads - 1;  // use double to avoid overflow
  return (x * sizeOfThisRange + sizeOfAllRanges) / sizeOfAllRanges;
}

void SubsetSuffixArray::sortRanges(std::vector<Range> *stacks,
				   indexT *bucketSizes,
				   const uchar *text,
				   size_t maxUnsortedInterval,
				   size_t numOfThreads) {
  std::vector<Range> &myStack = stacks[0];

  while( !myStack.empty() ){
#ifdef HAS_CXX_THREADS
    size_t numOfChunks = std::min(numOfThreads, myStack.size());
    if (numOfChunks > 1) {
      size_t totalSize = rangeSizeSum(myStack);
      size_t numOfNewThreads = numOfChunks - 1;
      std::vector<std::thread> threads(numOfNewThreads);

      for (size_t i = 0; i < numOfNewThreads; ++i) {
	size_t thisSize = nextRangeSize(myStack);
	size_t t = numOfThreadsForOneRange(numOfThreads, thisSize, totalSize);
	size_t maxThreads = numOfThreads - (numOfNewThreads - i);
	size_t thisThreads = std::min(t, maxThreads);
	numOfThreads -= thisThreads;

	do {
	  totalSize -= nextRangeSize(myStack);
	  stacks[numOfThreads].push_back(myStack.back());
	  myStack.pop_back();
	  thisSize += nextRangeSize(myStack);
	} while (myStack.size() > numOfThreads &&
		 thisSize <= totalSize / numOfThreads);
	// We want:
	// max(r | sizeSum(r) <= (totalSize - sizeSum(r-1)) / newNumOfThreads)

	threads[i] = std::thread(&SubsetSuffixArray::sortRanges, this,
				 stacks + numOfThreads,
				 bucketSizes + numOfThreads * numOfBuckets,
				 text, maxUnsortedInterval, thisThreads);
      }
      sortRanges(stacks, bucketSizes, text, maxUnsortedInterval, numOfThreads);
      for (size_t i = 0; i < numOfNewThreads; ++i) {
	threads[i].join();
      }
      return;
    }
#endif

    indexT* beg = myStack.back().beg;
    indexT* end = myStack.back().end;
    indexT depth = myStack.back().depth;
    myStack.pop_back();

    size_t interval = end - beg;
    const indexT minLength = 1;
    if( interval <= maxUnsortedInterval && depth >= minLength ) continue;

    const uchar* textBase = text + depth;
    const uchar* subsetMap = seed.subsetMap(depth);

    if( childTable.v.empty() && kiddyTable.v.empty() && chibiTable.v.empty() ){
      if( interval < 10 ){  // ???
	insertionSort( textBase, seed, beg, end, subsetMap );
	continue;
      }
    }else{
      if( interval == 2 ){
	sort2( textBase, beg, subsetMap );
	continue;
      }
    }

    unsigned subsetCount = seed.subsetCount(depth);

    ++depth;

    switch( subsetCount ){
    case 1:  radixSort1(myStack, textBase, subsetMap, beg, end, depth); break;
    case 2:  radixSort2(myStack, textBase, subsetMap, beg, end, depth); break;
    case 3:  radixSort3(myStack, textBase, subsetMap, beg, end, depth); break;
    case 4:  radixSort4(myStack, textBase, subsetMap, beg, end, depth); break;
    default: radixSortN(myStack, textBase, subsetMap, beg, end, depth,
			subsetCount, bucketSizes);
    }
  }
}

void SubsetSuffixArray::sortIndex( const uchar* text,
				   size_t maxUnsortedInterval,
				   int childTableType,
				   size_t numOfThreads ){
  if( childTableType == 1 ) chibiTable.v.assign( suffixArray.v.size(), -1 );
  if( childTableType == 2 ) kiddyTable.v.assign( suffixArray.v.size(), -1 );
  if( childTableType == 3 ) childTable.v.assign( suffixArray.v.size(), 0 );

  std::vector< std::vector<Range> > stacks(numOfThreads);
  std::vector<indexT> bucketSizes(numOfThreads * numOfBuckets);

  pushRange( stacks[0], &suffixArray.v.front(), &suffixArray.v.back() + 1, 0 );
  setChildReverse( &suffixArray.v.back() + 1, &suffixArray.v.front() );

  sortRanges(&stacks[0], &bucketSizes[0], text, maxUnsortedInterval,
	     numOfThreads);
}
