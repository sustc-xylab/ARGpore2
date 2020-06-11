CXXFLAGS = -O3 -Wall -Wextra -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -pedantic -g -std=c++11 -pthread -DHAS_CXX_THREADS
# -Wconversion
# -fomit-frame-pointer ?

CFLAGS = -Wall -O2

alpObj = alp/sls_alignment_evaluer.o alp/sls_pvalues.o		\
alp/sls_alp_sim.o alp/sls_alp_regression.o alp/sls_alp_data.o	\
alp/sls_alp.o alp/sls_basic.o alp/njn_localmaxstatmatrix.o	\
alp/njn_localmaxstat.o alp/njn_localmaxstatutil.o		\
alp/njn_dynprogprob.o alp/njn_dynprogprobproto.o		\
alp/njn_dynprogproblim.o alp/njn_ioutil.o alp/njn_random.o	\
alp/sls_falp_alignment_evaluer.o alp/sls_fsa1_pvalues.o		\
alp/sls_fsa1_utils.o alp/sls_fsa1.o alp/sls_fsa1_parameters.o

indexObj0 = Alphabet.o CyclicSubsetSeed.o LambdaCalculator.o	\
ScoreMatrix.o SubsetMinimizerFinder.o TantanMasker.o fileMap.o	\
cbrc_linalg.o mcf_substitution_matrix_stats.o tantan.o		\
LastdbArguments.o

indexObj4 = MultiSequence.o MultiSequenceQual.o SubsetSuffixArray.o	\
SubsetSuffixArraySort.o lastdb.o

alignObj0 = Alphabet.o CyclicSubsetSeed.o LambdaCalculator.o		\
ScoreMatrix.o SubsetMinimizerFinder.o TantanMasker.o fileMap.o		\
tantan.o LastalArguments.o GappedXdropAligner.o				\
GappedXdropAlignerPssm.o GappedXdropAligner2qual.o			\
GappedXdropAligner3frame.o mcf_gap_costs.o GeneticCode.o	\
GreedyXdropAligner.o LastEvaluer.o OneQualityScoreMatrix.o		\
QualityPssmMaker.o TwoQualityScoreMatrix.o gaplessXdrop.o		\
gaplessPssmXdrop.o gaplessTwoQualityXdrop.o cbrc_linalg.o		\
mcf_substitution_matrix_stats.o $(alpObj)

alignObj4 = Alignment.o AlignmentPot.o AlignmentWrite.o Centroid.o	\
DiagonalTable.o MultiSequence.o MultiSequenceQual.o SegmentPair.o	\
SegmentPairPot.o SubsetSuffixArray.o SubsetSuffixArraySearch.o		\
lastal.o

splitObj0 = Alphabet.o LambdaCalculator.o fileMap.o cbrc_linalg.o	\
mcf_substitution_matrix_stats.o split/cbrc_unsplit_alignment.o		\
split/last-split-main.o

splitObj4 = MultiSequence.o split/cbrc_split_aligner.o	\
split/last-split.o

PPOBJ = last-pair-probs.o last-pair-probs-main.o

MBOBJ = last-merge-batches.o

ALL = lastdb lastal last-split last-merge-batches last-pair-probs	\
lastdb8 lastal8 last-split8

indexObj8 = $(indexObj4:.o=.o8)
alignObj8 = $(alignObj4:.o=.o8)
splitObj8 = $(splitObj4:.o=.o8)

all: $(ALL)

indexAllObj4 = $(indexObj0) $(indexObj4)
lastdb: $(indexAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexAllObj4) -lz

indexAllObj8 = $(indexObj0) $(indexObj8)
lastdb8: $(indexAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexAllObj8) -lz

alignAllObj4 = $(alignObj0) $(alignObj4)
lastal: $(alignAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignAllObj4) -lz

alignAllObj8 = $(alignObj0) $(alignObj8)
lastal8: $(alignAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignAllObj8) -lz

splitAllObj4 = $(splitObj0) $(splitObj4)
last-split: $(splitAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitAllObj4)

splitAllObj8 = $(splitObj0) $(splitObj8)
last-split8: $(splitAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitAllObj8)

last-pair-probs: $(PPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ) -lz

last-merge-batches: $(MBOBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MBOBJ)

.SUFFIXES:
.SUFFIXES: .o .c .cc .cpp .o8

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

.cc.o:
	$(CXX) -DLAST_INT_TYPE=unsigned $(CPPFLAGS) $(CXXFLAGS) -I. -c -o $@ $<

.cc.o8:
	$(CXX) -DLAST_INT_TYPE=size_t $(CPPFLAGS) $(CXXFLAGS) -I. -c -o $@ $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(ALL) *.o* */*.o*

CyclicSubsetSeedData.hh: ../data/*.seed
	../build/seed-inc.sh ../data/*.seed > $@

ScoreMatrixData.hh: ../data/*.mat
	../build/mat-inc.sh ../data/*.mat > $@

VERSION = \"`hg id -n`\"
VERSION2 = \"`grep latesttagdistance ../.hg_archival.txt | cut -d' ' -f2`\"
VERSION3 = \"UNKNOWN\"

version.hh: FORCE
	if test -e ../.hg ; \
	then echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@ ; \
	elif test -e ../.hg_archival.txt ; \
	then echo $(VERSION2) | cmp -s $@ - || echo $(VERSION2) > $@ ; \
	else test -e $@ || echo $(VERSION3) > $@ ; \
	fi

FORCE:

fix8 = sed 's/.*\.o/& &8/'

depend:
	sed '/[m][v]/q' makefile > m
	$(CXX) -MM *.cc | $(fix8) >> m
	$(CC) -MM *.c >> m
	$(CXX) -MM alp/*.cpp | sed 's|.*:|alp/&|' >> m
	$(CXX) -MM -I. split/*.cc | sed 's|.*:|split/&|' | $(fix8) >> m
	mv m makefile
Alignment.o Alignment.o8: Alignment.cc Alignment.hh ScoreMatrixRow.hh SegmentPair.hh \
 mcf_gap_costs.hh Alphabet.hh Centroid.hh GappedXdropAligner.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh GeneticCode.hh \
 GreedyXdropAligner.hh TwoQualityScoreMatrix.hh
AlignmentPot.o AlignmentPot.o8: AlignmentPot.cc AlignmentPot.hh Alignment.hh \
 ScoreMatrixRow.hh SegmentPair.hh mcf_gap_costs.hh
AlignmentWrite.o AlignmentWrite.o8: AlignmentWrite.cc Alignment.hh ScoreMatrixRow.hh \
 SegmentPair.hh mcf_gap_costs.hh GeneticCode.hh LastEvaluer.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 MultiSequence.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 Alphabet.hh
Alphabet.o Alphabet.o8: Alphabet.cc Alphabet.hh
Centroid.o Centroid.o8: Centroid.cc Centroid.hh GappedXdropAligner.hh \
 ScoreMatrixRow.hh mcf_gap_costs.hh SegmentPair.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 GappedXdropAlignerInl.hh
CyclicSubsetSeed.o CyclicSubsetSeed.o8: CyclicSubsetSeed.cc CyclicSubsetSeed.hh \
 CyclicSubsetSeedData.hh zio.hh mcf_zstream.hh stringify.hh
DiagonalTable.o DiagonalTable.o8: DiagonalTable.cc DiagonalTable.hh
GappedXdropAligner.o GappedXdropAligner.o8: GappedXdropAligner.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner2qual.o GappedXdropAligner2qual.o8: GappedXdropAligner2qual.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh \
 TwoQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh
GappedXdropAligner3frame.o GappedXdropAligner3frame.o8: GappedXdropAligner3frame.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner3framePssm.o GappedXdropAligner3framePssm.o8: GappedXdropAligner3framePssm.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerPssm.o GappedXdropAlignerPssm.o8: GappedXdropAlignerPssm.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GeneticCode.o GeneticCode.o8: GeneticCode.cc GeneticCode.hh Alphabet.hh
GreedyXdropAligner.o GreedyXdropAligner.o8: GreedyXdropAligner.cc GreedyXdropAligner.hh \
 ScoreMatrixRow.hh
LambdaCalculator.o LambdaCalculator.o8: LambdaCalculator.cc LambdaCalculator.hh \
 cbrc_linalg.hh
LastEvaluer.o LastEvaluer.o8: LastEvaluer.cc LastEvaluer.hh ScoreMatrixRow.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 GeneticCode.hh
LastalArguments.o LastalArguments.o8: LastalArguments.cc LastalArguments.hh \
 SequenceFormat.hh stringify.hh getoptUtil.hh version.hh
LastdbArguments.o LastdbArguments.o8: LastdbArguments.cc LastdbArguments.hh \
 SequenceFormat.hh stringify.hh getoptUtil.hh version.hh
MultiSequence.o MultiSequence.o8: MultiSequence.cc MultiSequence.hh ScoreMatrixRow.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh io.hh
MultiSequenceQual.o MultiSequenceQual.o8: MultiSequenceQual.cc MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
OneQualityScoreMatrix.o OneQualityScoreMatrix.o8: OneQualityScoreMatrix.cc \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
QualityPssmMaker.o QualityPssmMaker.o8: QualityPssmMaker.cc QualityPssmMaker.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
ScoreMatrix.o ScoreMatrix.o8: ScoreMatrix.cc ScoreMatrix.hh ScoreMatrixData.hh \
 qualityScoreUtil.hh stringify.hh zio.hh mcf_zstream.hh
SegmentPair.o SegmentPair.o8: SegmentPair.cc SegmentPair.hh
SegmentPairPot.o SegmentPairPot.o8: SegmentPairPot.cc SegmentPairPot.hh SegmentPair.hh
SubsetMinimizerFinder.o SubsetMinimizerFinder.o8: SubsetMinimizerFinder.cc \
 SubsetMinimizerFinder.hh CyclicSubsetSeed.hh
SubsetSuffixArray.o SubsetSuffixArray.o8: SubsetSuffixArray.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 SubsetMinimizerFinder.hh io.hh
SubsetSuffixArraySearch.o SubsetSuffixArraySearch.o8: SubsetSuffixArraySearch.cc \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh
SubsetSuffixArraySort.o SubsetSuffixArraySort.o8: SubsetSuffixArraySort.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
TantanMasker.o TantanMasker.o8: TantanMasker.cc TantanMasker.hh ScoreMatrixRow.hh \
 tantan.hh ScoreMatrix.hh mcf_substitution_matrix_stats.hh
TwoQualityScoreMatrix.o TwoQualityScoreMatrix.o8: TwoQualityScoreMatrix.cc \
 TwoQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
cbrc_linalg.o cbrc_linalg.o8: cbrc_linalg.cc cbrc_linalg.hh
fileMap.o fileMap.o8: fileMap.cc fileMap.hh stringify.hh
gaplessPssmXdrop.o gaplessPssmXdrop.o8: gaplessPssmXdrop.cc gaplessPssmXdrop.hh \
 ScoreMatrixRow.hh
gaplessTwoQualityXdrop.o gaplessTwoQualityXdrop.o8: gaplessTwoQualityXdrop.cc \
 gaplessTwoQualityXdrop.hh TwoQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh ScoreMatrixRow.hh
gaplessXdrop.o gaplessXdrop.o8: gaplessXdrop.cc gaplessXdrop.hh ScoreMatrixRow.hh
last-pair-probs-main.o last-pair-probs-main.o8: last-pair-probs-main.cc last-pair-probs.hh \
 stringify.hh version.hh
last-pair-probs.o last-pair-probs.o8: last-pair-probs.cc last-pair-probs.hh zio.hh \
 mcf_zstream.hh stringify.hh
lastal.o lastal.o8: lastal.cc last.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 SequenceFormat.hh qualityScoreUtil.hh LastalArguments.hh \
 QualityPssmMaker.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh TwoQualityScoreMatrix.hh LastEvaluer.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 GeneticCode.hh SubsetMinimizerFinder.hh SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh Centroid.hh GappedXdropAligner.hh mcf_gap_costs.hh \
 SegmentPair.hh AlignmentPot.hh Alignment.hh SegmentPairPot.hh \
 ScoreMatrix.hh TantanMasker.hh tantan.hh DiagonalTable.hh \
 GreedyXdropAligner.hh gaplessXdrop.hh gaplessPssmXdrop.hh \
 gaplessTwoQualityXdrop.hh zio.hh mcf_zstream.hh threadUtil.hh version.hh
lastdb.o lastdb.o8: lastdb.cc last.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 SequenceFormat.hh qualityScoreUtil.hh LastdbArguments.hh \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh TantanMasker.hh tantan.hh \
 zio.hh mcf_zstream.hh threadUtil.hh version.hh
mcf_gap_costs.o mcf_gap_costs.o8: mcf_gap_costs.cc mcf_gap_costs.hh
mcf_substitution_matrix_stats.o mcf_substitution_matrix_stats.o8: mcf_substitution_matrix_stats.cc \
 mcf_substitution_matrix_stats.hh LambdaCalculator.hh cbrc_linalg.hh
tantan.o tantan.o8: tantan.cc tantan.hh
last-merge-batches.o: last-merge-batches.c version.hh
alp/njn_dynprogprob.o: alp/njn_dynprogprob.cpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogproblim.o: alp/njn_dynprogproblim.cpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogprobproto.o: alp/njn_dynprogprobproto.cpp \
 alp/njn_dynprogprobproto.hpp
alp/njn_ioutil.o: alp/njn_ioutil.cpp alp/njn_ioutil.hpp
alp/njn_localmaxstat.o: alp/njn_localmaxstat.cpp alp/sls_basic.hpp \
 alp/njn_localmaxstat.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_function.hpp alp/njn_doubletype.hpp \
 alp/njn_integer.hpp alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp \
 alp/njn_approx.hpp alp/njn_vector.hpp
alp/njn_localmaxstatmatrix.o: alp/njn_localmaxstatmatrix.cpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_memutil.hpp
alp/njn_localmaxstatutil.o: alp/njn_localmaxstatutil.cpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_integer.hpp alp/njn_memutil.hpp \
 alp/njn_root.hpp alp/njn_function.hpp alp/sls_basic.hpp
alp/njn_random.o: alp/njn_random.cpp alp/njn_random.hpp
alp/sls_alignment_evaluer.o: alp/sls_alignment_evaluer.cpp \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp.hpp alp/sls_alp_data.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_sim.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_alp.o: alp/sls_alp.cpp alp/sls_alp.hpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_data.o: alp/sls_alp_data.cpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_regression.o: alp/sls_alp_regression.cpp \
 alp/sls_alp_regression.hpp alp/sls_basic.hpp
alp/sls_alp_sim.o: alp/sls_alp_sim.cpp alp/sls_alp_sim.hpp \
 alp/sls_alp_data.hpp alp/sls_basic.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp.hpp
alp/sls_basic.o: alp/sls_basic.cpp alp/sls_basic.hpp
alp/sls_falp_alignment_evaluer.o: alp/sls_falp_alignment_evaluer.cpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/sls_fsa1.hpp
alp/sls_fsa1.o: alp/sls_fsa1.cpp alp/sls_fsa1.hpp alp/sls_alp_regression.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_pvalues.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_fsa1_parameters.o: alp/sls_fsa1_parameters.cpp \
 alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp alp/sls_basic.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/sls_fsa1_pvalues.hpp
alp/sls_fsa1_pvalues.o: alp/sls_fsa1_pvalues.cpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
alp/sls_fsa1_utils.o: alp/sls_fsa1_utils.cpp alp/sls_fsa1_utils.hpp \
 alp/sls_basic.hpp alp/njn_random.hpp alp/njn_uniform.hpp
alp/sls_pvalues.o: alp/sls_pvalues.cpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp_data.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
split/cbrc_split_aligner.o split/cbrc_split_aligner.o8: split/cbrc_split_aligner.cc \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 mcf_substitution_matrix_stats.hh
split/cbrc_unsplit_alignment.o split/cbrc_unsplit_alignment.o8: split/cbrc_unsplit_alignment.cc \
 split/cbrc_unsplit_alignment.hh
split/last-split-main.o split/last-split-main.o8: split/last-split-main.cc split/last-split.hh \
 stringify.hh version.hh
split/last-split.o split/last-split.o8: split/last-split.cc split/last-split.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
