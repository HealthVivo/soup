SAMTOOLS=samtools-0.1.19
CFLAGS = -I ${SAMTOOLS} -L${SAMTOOLS} -lbam -lz -O3

bamtofastq: bamtofastq.cpp
	g++ -o bamtofastq bamtofastq.cpp ${CFLAGS}


