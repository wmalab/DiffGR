# DiffGR

DiffGR is a novel statistical method for detecting differential genomic regions at TAD level between two Hi-C contact maps. Briefly, DiffGR utilizes the stratum-adjusted correlation coefficient (SCC), which can effectively eliminate the genomic-distance effect in Hi-C data, to measure the similarity of local genomic regions between two contact matrices, and then applies a nonparametric permutation test on those SCC values to detect genomic regions. 


# Installation

The source code can be performed under R language version 4.0.2 with the installation of packages HiCcompare, HiCseg and R.utils.



# Input

```

dat1,dat2             numeric.  N*N raw HiC contact maps, which would firstly be preprocessed with 2D
                      mean filter smoothing and KR normalization in DiffGR function for the later use

tad1,tad2             numeric. A vector of TAD boundaries of contact maps.If the input is NA, the program
                      will automatically detect the TADs by HiCseg

res                   numeric. The resolution of HiC contact maps, eg:100kb will input 100,000

smooth.size           numeric. The size controlling the smoothing level (The size varies across different
                      resolution and is guided by Hicrep paper). In this paper, we obtained the smoothing
                      size with 11, 5 and 3 on real data analysis for the resolution of 25Kb, 50Kb and 
                      100Kb respectively, and set the smoothing size with 0 in simulation.
                      
N.perm                numeric. The number of iterations in permutation test

cutoff.default        logical. Whether set the SCC cutoff (meaningful SCC between the two Hi-C datasets
                      that mustbe reached in order to call a differential TAD truly significant) with
                      self-defined value(True) or with automatic computed value (False)
                      
speedup.option        logical. Calculation with or without speed-up algorithm (True/FALSE)

alpha                 numeric. Significant level of differential region testing 

```

# Output

return a list that contains the tad result and genomic region result


The tad result table contains the following elements:

```

tad.start              the start location for the starting bin of TAD

tad.end                the start location for the end bin of TAD

scc                    the SCC value of corresponding domain

pvalue                 the pvalue of differential testing on corresponding domain

pvalue.adj             the adjusted pvalue of differential testing on corresponding domain 
                       (adjusted by Benjamin-Hochberg)
                       
```

The genomic result table contains the following elements:

```

genom.start            tthe start location for the starting bin of genomic region

genom.end              the start location for the end bin of genomic region

condition.type         the type if candidate genomic region belonging to. 
                       1:single-TAD, 2: Hierachical-TAD, 3: Alternating-TAD

detect.result          the differential testing result for corresponding genomic region. 
                       1:Differential 0:Non-differential 
                       
```

# Sample Data

The raw HiC contact maps getting from chr10 of GM12878 and HMEC with resolution=50kb were untilized as sample data. An example of the usage of DiffGR with/without TAD inputs is shown below:

```
dat1 <- readRDS("path/dat.GM12878.chr10.rds")
dat2 <- readRDS("path/dat.K562.chr10.rds")
tad1 <- read.table("path/tad.GM12878.chr10.txt")
tad1 <- tad1$x
tad2 <- read.table("path/tad.K562.chr10.txt") 
tad2 <- tad2$x

#with TAD inputs
result <- DiffGR(dat1=dat1,dat2=dat2,tad1=tad1,tad2=tad2,smooth.size=5,res=50000)

#without TAD inputs
result <- DiffGR(dat1=dat1,dat2=dat2,smooth.size=5,res=50000)
```



