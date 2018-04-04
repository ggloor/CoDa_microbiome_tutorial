This directory contains the readings, materials, and examples for a workshop originally offered at the Exploring Human Host-Microbiome Interactions in Health and Disease 2016 conference.

The workshop was subsequently modified for the NGS'18 conference in Barcelona. For this, we examine an RNA-seq dataset, outlier detection was added, and the use of the propr R package was added (Ionas Erb).

The primary focus will be on explaining why transcriptome data are compositional, limits and caveats on the analysis of compositional data, and finally proper methods to analyze such data. Practical methods to deal with the sparse nature of the data will be presented and demonstrated in practice.

The workshop will focus on hands-on application of three main methods.

Exploratory data analysis using the compositional biplot which shows the variance in both the samples and the genes on one plot. Outlier detection and removal

Proper correlation and association analysis using the latest methods including the rho-metric inside propr and other tools for correlation.

Examination of differential abundance under a compositional approach using the ALDEx2 tool.

Workshop participants should read the following:

multi_comp.pdf, and coda_seq.pdf: brief introductions to the problem from different perspectives
make_interpret_biplot: a simple biplot with random data
first_biplot.pdf: exploratory data analysis using compositional biplots and outlier detection and removal of the example dataset
first_comparison.pdf: finding maximally different genes using ALDEx2
first_association.pdf: finding maximally associated genes using propr
