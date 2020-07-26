DiffGR: Differential genomic region detection function




@dat1,@dat2 numeric. N*N HiC contact maps which have been preprocessed (2D mean filter smoothing and KR normalization)

@tad1 @tad2 numeric. TAD boundaries of dat1 and dat2 (calling by HiCseg)

@N.perm numeric. The number of iterations in permutation test

@speedup.option logical. (True/FALSE) Calculation with or without speed-up algorithm

@alpha logical. Significant level of differential region testing 




@return a list that contains the tad result and genomic region result


The tad result table contains the following elements:

tad.start: the starting locus of TAD

tad.end: the end locus of TAD

scc: the SCC value of corresponding domain

pvalue: the pvalue of differential testing on corresponding domain

pvalue.adj: the adjusted pvalue of differential testing on corresponding domain (adjusted by Benjamin-Hochberg)


The genomic result table contains the following elements:

genom.start: the starting locus of genomic region

genom.end: the end locus of genomic region

condition.type: the type if candidate genomic region belonging to. 1:single-TAD, 2: Hierachical-TAD, 3: Alternating-TAD

detect.result: the differential testing result for corresponding genomic region. 1:Differential 0:Non-differential 

