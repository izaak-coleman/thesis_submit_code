
SMuFin: main.cpp alphabet.o amino.o array8.o array.o direct.o genarray.o genhash.o gmem.o gs2k.o hash.o mathutil.o mystring.o parse.o seqgen.o stack.o stree.o trace.o util.o blind2.o deep2.o ds.o globals.o helped.o shallow.o getway_s asm bwa.a samtools.a 
	mpic++ -o $@ main.cpp alphabet.o amino.o array8.o array.o direct.o genarray.o genhash.o gmem.o gs2k.o hash.o mathutil.o mystring.o parse.o seqgen.o stack.o stree.o trace.o util.o blind2.o deep2.o ds.o globals.o helped.o shallow.o getway_s asm bwa.a samtools.a -lz -lmpi -lgomp -std=c++11

bwtindex.o: bwa/bwtindex.c
	gcc -fopenmp -c bwa/bwtindex.c

utils.o: bwa/utils.c
	gcc -fopenmp -c bwa/utils.c

is.o: bwa/is.c
	gcc -fopenmp -c bwa/is.c

bwt.o: bwa/bwt.c
	gcc -fopenmp -c bwa/bwt.c

bntseq.o: bwa/bntseq.c
	gcc -fopenmp -c bwa/bntseq.c

bwt_gen.o: bwa/bwt_gen.c
	gcc -fopenmp -c bwa/bwt_gen.c

QSufSort.o: bwa/QSufSort.c
	gcc -fopenmp -c bwa/QSufSort.c

bwtaln.o: bwa/bwtaln.c
	gcc -fopenmp -c bwa/bwtaln.c

bwtgap.o: bwa/bwtgap.c
	gcc -fopenmp -c bwa/bwtgap.c

bwaseqio.o: bwa/bwaseqio.c
	gcc -fopenmp -c bwa/bwaseqio.c

bwa.o: bwa/bwa.c
	gcc -fopenmp -c bwa/bwa.c

ksw.o: bwa/ksw.c
	gcc -fopenmp -c bwa/ksw.c

bwase.o: bwa/bwase.c
	gcc -fopenmp -c bwa/bwase.c

bwape.o: bwa/bwape.c
	gcc -fopenmp -c bwa/bwape.c

bwamem.o: bwa/bwamem.c
	gcc -fopenmp -c bwa/bwamem.c

bwamem_pair.o: bwa/bwamem_pair.c
	gcc -fopenmp -c bwa/bwamem_pair.c

fastmap.o: bwa/fastmap.c
	gcc -fopenmp -c bwa/fastmap.c

kopen.o: bwa/kopen.c
	gcc -fopenmp -c bwa/kopen.c

kstring.o: bwa/kstring.c
	gcc -fopenmp -c bwa/kstring.c

sam_view.o: samtools-0.1.19/sam_view.c
	gcc -fopenmp -c samtools-0.1.19/sam_view.c

sam_header.o: samtools-0.1.19/sam_header.c
	gcc -fopenmp -c samtools-0.1.19/sam_header.c

sam.o: samtools-0.1.19/sam.c
	gcc -fopenmp -c samtools-0.1.19/sam.c

knetfile.o: samtools-0.1.19/knetfile.c
	gcc -fopenmp -c samtools-0.1.19/knetfile.c

bam.o: samtools-0.1.19/bam.c
	gcc -fopenmp -c samtools-0.1.19/bam.c

bam_sort.o: samtools-0.1.19/bam_sort.c
	gcc -fopenmp -c samtools-0.1.19/bam_sort.c

bam_index.o: samtools-0.1.19/bam_index.c
	gcc -fopenmp -c samtools-0.1.19/bam_index.c

bgzf.o: samtools-0.1.19/bgzf.c
	gcc -fopenmp -c samtools-0.1.19/bgzf.c

bam_aux.o: samtools-0.1.19/bam_aux.c
	gcc -fopenmp -c samtools-0.1.19/bam_aux.c

bam_import.o: samtools-0.1.19/bam_import.c
	gcc -fopenmp -c samtools-0.1.19/bam_import.c

bam_pileup.o: samtools-0.1.19/bam_pileup.c
	gcc -fopenmp -c samtools-0.1.19/bam_pileup.c

faidx.o: samtools-0.1.19/faidx.c
	gcc -fopenmp -c samtools-0.1.19/faidx.c

razf.o: samtools-0.1.19/razf.c
	gcc -fopenmp -c samtools-0.1.19/razf.c

bedidx.o: samtools-0.1.19/bedidx.c
	gcc -fopenmp -c samtools-0.1.19/bedidx.c

alphabet.o: alphabet.c
	gcc -c alphabet.c

amino.o: amino.c
	gcc -c amino.c

array8.o: array8.c
	gcc -c array8.c

array.o: array.c
	gcc -c array.c

direct.o: direct.c
	gcc -c direct.c

genarray.o: genarray.c
	gcc -c genarray.c

genhash.o: genhash.c
	gcc -c genhash.c

gmem.o: gmem.c
	gcc -c gmem.c

gs2k.o: gs2k.c
	gcc -c gs2k.c

hash.o: hash.c
	gcc -c hash.c

mathutil.o: mathutil.c
	gcc -c mathutil.c

mystring.o: mystring.c
	gcc -c mystring.c

parse.o: parse.c
	gcc -c parse.c

seqgen.o: seqgen.c
	gcc -c seqgen.c

stack.o: stack.c
	gcc -c stack.c

stree.o: stree.c
	gcc -c stree.c

trace.o: trace.c
	gcc -c trace.c

util.o: util.c
	gcc -c util.c

blind2.o: blind2.c
	gcc -c blind2.c

deep2.o: deep2.c
	gcc -c deep2.c

ds.o: ds.c
	gcc -c ds.c

globals.o: globals.c
	gcc -c globals.c

helped.o: helped.c
	gcc -c helped.c

shallow.o: shallow.c
	gcc -c shallow.c

bwa.a: bwtindex.o utils.o is.o bwt.o bntseq.o bwt_gen.o QSufSort.o bwtaln.o bwtgap.o bwaseqio.o bwa.o ksw.o bwase.o bwape.o kstring.o bwamem.o bwamem_pair.o fastmap.o kopen.o
	ar rcs bwa.a bwtindex.o utils.o is.o bwt.o bntseq.o bwt_gen.o QSufSort.o bwtaln.o bwtgap.o bwaseqio.o bwa.o ksw.o bwase.o bwape.o kstring.o bwamem.o bwamem_pair.o fastmap.o kopen.o kstring.o

samtools.a: sam_view.o sam_header.o sam.o knetfile.o bam.o bam_sort.o bam_index.o bgzf.o bam_aux.o bam_import.o bam_pileup.o faidx.o razf.o kstring.o bedidx.o
	ar rcs samtools.a sam_view.o sam_header.o sam.o knetfile.o bam.o bam_sort.o bam_index.o bgzf.o bam_aux.o bam_import.o bam_pileup.o faidx.o razf.o kstring.o bedidx.o

clean:
	rm -f *.o *.a SMuFin
