CC=/usr/bin/gcc -O0 -g

SOBJ = main.o sw.o utils.o hash.o hash_func.o str.o gap_closer.o kmer.o contig.o ont.o rseq.o bio.o gc_graph.o digraph.o cigar.o ctg_graph.o lfr.o

ont: $(SOBJ)
	$(CC) -o gc $(SOBJ) -lm -lz -lpthread

.PHONY: clean
clean:
	rm -f *.o gc
