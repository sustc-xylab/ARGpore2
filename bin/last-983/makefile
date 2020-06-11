CXXFLAGS = -O3 -std=c++11 -pthread -DHAS_CXX_THREADS
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

progs = src/lastdb src/lastal src/last-split src/last-merge-batches	\
src/last-pair-probs src/lastdb8 src/lastal8 src/last-split8

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp $(progs) scripts/* $(bindir)

clean:
	@cd src && $(MAKE) clean

html:
	@cd doc && $(MAKE)

distdir = last-`hg id -n`

RSYNCFLAGS = -aC --exclude '*8' --exclude 'last??' --exclude last-split --exclude last-merge-batches --exclude last-pair-probs

dist: log html
	@cd src && $(MAKE) version.hh CyclicSubsetSeedData.hh ScoreMatrixData.hh
	rsync $(RSYNCFLAGS) build doc examples makefile scripts src data *.txt $(distdir)
	zip -qrm $(distdir) $(distdir)

log:
	hg log --style changelog > ChangeLog.txt
