# DiffGR

DiffGR is a novel statistical method for detecting differential genomic regions at TAD level between two Hi-C contact maps. Briefly, DiffGR utilizes the stratum-adjusted correlation coefficient (SCC), which can effectively eliminate the genomic-distance effect in Hi-C data, to measure the similarity of local genomic regions between two contact matrices, and then applies a nonparametric permutation test on those SCC values to detect genomic regions. 


# Installation

The source code can be performed under R language version 4.0.2 with the installation of packages HiCcompare, hicrep，HiCseg and R.utils.



# Input

```

@dat1,@dat2 numeric. N*N HiC contact maps which have been preprocessed with 2D mean filter smoothing and KR normalization

@res numeric. The resolution of HiC contact maps, eg:100kb will input 100,000

@N.perm numeric. The number of iterations in permutation test

@speedup.option logical. (True/FALSE) Calculation with or without speed-up algorithm

@alpha logical. Significant level of differential region testing 

```

# Output

return a list that contains the tad result and genomic region result


The tad result table contains the following elements:

```

tad.start: the starting locus of TAD

tad.end: the end locus of TAD

scc: the SCC value of corresponding domain

pvalue: the pvalue of differential testing on corresponding domain

pvalue.adj: the adjusted pvalue of differential testing on corresponding domain (adjusted by Benjamin-Hochberg)
```

The genomic result table contains the following elements:

```

genom.start: the starting locus of genomic region

genom.end: the end locus of genomic region

condition.type: the type if candidate genomic region belonging to. 1:single-TAD, 2: Hierachical-TAD, 3: Alternating-TAD

detect.result: the differential testing result for corresponding genomic region. 1:Differential 0:Non-differential 
```

# Sample Data

The data getting from chr10 of GM12878 and HMEC with resolution=50kb were untilized as sample data. In the sample data file, dat1 and dat2 represent their corresponding HiC contacts with 2D mean filter smoothing and KR normalization; tad1 and tad2 denote their TAD boundaries which were detected by HiCseg. 

To run the sampe data, 

```
dat1 <- readRDS("path/dat1.rds")
dat2 <- readRDS("path/dat2.rds")
tad1 <- read.table("path/tad1.txt")
tad2 <- read.table("path/tad2.txt")
result <- DiffGR(dat1,tad1$x,dat2,tad2$x,res=50000)

```



