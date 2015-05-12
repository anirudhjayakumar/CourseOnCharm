CHARMC = /home/ajayaku2/charm/netlrts-linux-x86_64/bin/charmc $(OPTS)
TESTOPTS = ++local

OBJS = particle.o

N = 100
K = 16

all: particle projection

particle: $(OBJS)
	$(CHARMC) -O3 -language charm++ -module liveViz -module CommonLBs -o particle $(OBJS)

projection: $(OBJS)
	$(CHARMC) -O3 -language charm++ -module liveViz -module CommonLBs -tracemode projections -lz -o particle.prj $(OBJS)
	
particle.decl.h: particle.ci
	$(CHARMC)  particle.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o particle charmrun particle.prj

cleanp:
	rm -f *.sts *.gz *.projrc *.topo *.out

particle.o: particle.C particle.decl.h
	$(CHARMC) -c -g particle.C

test: all
	./charmrun particle $(N) $(K) $(TESTOPTS) +p4
