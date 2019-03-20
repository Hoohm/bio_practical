library(tidyr)
library(DESeq2)
library(gtools)
library(ggplot2)
library(dplyr)

# There is a tutorial that goes deep into Deseq2 here: bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#First we set our working directory. We set it to the folder where our data is
setwd('')

#LOADING THE DATA
# Here we load the two files. The count matrix and the design file.
#Exercise: Find the right arguments to load the data to be a dataframe with samples in the colnames and genes in the rownames like below:

# L1-Seq34_S1 L2-Seq34_S2 L3-Seq34_S3 L4-Seq34_S4 L5-Seq34_S5 L6-Seq34_S6 L7-Seq34_S7 L8-Seq34_S8 L9-Seq34_S9 L10-Seq34_S10
# 0610005C13Rik           1           1           2           1           6           1           2           2           2             1
# 0610006L08Rik           0           0           0           0           0           0           0           0           0             0
# 0610009B22Rik         209         185         210         246         259         318         240         360         267           303
# 0610009E02Rik           1           2           2           3           3           5           2           0           0             7
# 0610009L18Rik          37          48          55          63          55          70          29          59          59            54
# 0610009O20Rik         142         132         162         194         165         214         146         219         172           169
expression = read.table()
#Exercise: Find the right arguments to load the data to be a dataframe with samples in the rownames and conditions in the colnames like below:
# Day Infection
# L1-Seq34_S1    D8       Arm
# L2-Seq34_S2    D8       Arm
# L3-Seq34_S3    D8       Arm
# L4-Seq34_S4    D8       Arm
# L5-Seq34_S5    D8       Arm
# L6-Seq34_S6    D8       c13
# L7-Seq34_S7    D8       c13
design = read.table()

# To be sure that we have the same order in both files we get a sorted vector of our samples
sorted_samples = mixedsort(rownames(design))

#We reorder the two dataframes using our vector.
expression = expression[,sorted_samples]
design = design[sorted_samples,]

#To simplify the design we combine both our contitions together.
#unite makes a combined version of our two columns
combined = unite(design, 'groups', c('Day','Infection'))

#To create a Deseq object, we need a design formula. This will tell the function what kind of differences we are looking for.
design_formula = as.formula('~ groups')

#Now we can create our deseq object using the count dataframe, design dataframe and the formula.
dds = DESeqDataSetFromMatrix(expression, colData = combined, design = design_formula)

#We filter lowly epxressed genes. This is mostly for performances purposes.
# Since we use a LOGICAL operator here, we will get ids back. Those ids correspond to the row number of the rows that fit the criteria we are using.
#Here we look for genes that are expressed  with a value lower or equal to the number of samples we have.

keep <- which(rowSums(counts(dds)) >= ncol(counts(dds)))
dds <- dds[keep,]

#Before we go to differential gene expression we want to check if our samples look fine.
#For this we use Principal Component Analysis (PCA).
#We need to normalize the data before we do this. Hence we use the Variance Stabilizing Transformation function here or vst.

#To have an idea of what the transformation does you can quickly make boxplots of the data before and after transformation

#Before
boxplot(counts(dds))
plot(density(counts(dds)))
#Do the transformation
vst = vst(dds)

#After
boxplot(assay(vst))
plot(density(assay(vst)))

#Notice the difference in the value range. This is because we are now using log values of the counts.

#To make a PCA we are not going to use all the genes. We are interested in genes that are highly variable accross the samples because they will give us the most information about our data.

#To decide which genes are more variant we can use different methods

#Mean absolute deviation
rm = rowMads(assay(vst))
#Interquantile range
riqr = rowIQRs(assay(vst))
#variance
rv = rowVars(assay(vst))
#All these methods are calculating how variable the expression of genes are across the samples.

#The rowmad selection is already done for you. Try to make the IQR and RV selection


################### MAD selection
#Before running this, look at the rm vector and the summary(rm) output. Can you tell what the summary function does?
plot(sort(rm))
boxplot(rm)

select <- which(rm > summary(rm)[5])

PCA <- prcomp(t(assay(vst)[select, ]), scale. = TRUE, center = TRUE)
plot_data = cbind(PCA$x, design, combined)

ggplot(plot_data, aes(x=PC1, y=PC2, group=groups, color=groups, label=rownames(plot_data))) + stat_ellipse(type='norm') + geom_point() + geom_text_repel(size=3) + ggtitle('PCA Dim1 and Dim2 using MAD for gene selection')
###################################



###################Exercise: IQR based selection

#########################

##################Exercise: Variance based selection

#####################

#PCA doesn't give us only 2 dimentions. can you create the same plot than before but with dimension 1 and 3 and another plot with dimension 2 and 3.?

######################Exercise: PCA plot 1-3

###################


######################Exercise: PCA plot 2-3

########################


########### Differential Expression Analysis
#### Now that we have seen that our data has no outliers we can proceed to the DGE.
dds = DESeq(dds)


res <- results(dds, contrast = c('groups', 'D8_Mix', 'D8_Arm'))


###Exercies: Find the differentialy expressed genes between Day 32 Arm and Day 32 c13

###############################






source('volcano.R')

