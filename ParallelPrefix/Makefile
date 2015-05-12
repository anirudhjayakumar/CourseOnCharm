OPTS = -g
CHARMC=/dcsdata/home/anirudh/code/charm/netlrts-linux-x86_64/bin/charmc $(OPTS)
#CHARMC=charmc $(OPTS)

all: ldb

ldb:   ldb.o
	$(CHARMC) ldb.o -o ldb -language charm++ 

ldb.o : ldb.cc ldb.def.h ldb.decl.h
	$(CHARMC)  -c ldb.cc

ldb.def.h ldb.decl.h : ldb.ci
	$(CHARMC) ldb.ci

test: all
	./charmrun ldb 8 100 1000 +p4 ++local

clean:
	rm -f conv-host *.o charmrun ldb
	rm -f *.def.h *.decl.h
	rm -f *.log.gz *.projrc *.topo *.sts

