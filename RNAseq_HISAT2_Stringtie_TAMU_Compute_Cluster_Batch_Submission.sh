#######################
## RNA-seq
#######################
for filename in /scratch/user/brian.a.mckinley/nroots/*.fastq.gz
do
cat > ${filename}.sh << EOF
#BSUB -J starch3
#BSUB -L /bin/bash
#BSUB -W 04:00
#BSUB -n 20
#BSUB -R "span[ptile=10]"
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000
#BSUB -o starch2_stdout.%J
#BSUB -e starch2_stderr.%J
module load HISAT2/2.0.5-intel-2015B-Python-2.7.10
module load SAMtools/1.3-intel-2015B
module load StringTie/1.3.3b-intel-2015B
hisat2 -p 20 -q -x /scratch/user/brian.a.mckinley/index/sbi_base -U ${filename} -S ${filename}.sam --known-splicesite-infile /scratch/user/brian.a.mckinley/referanceV3/Sbicolorsplicesites.txt;
samtools view -@ 20 -bS ${filename}.sam > ${filename}.bam;
samtools sort -@ 20 ${filename}.bam -o ${filename}.sorted.bam;
stringtie ${filename}.sorted.bam -p 20 -G /scratch/user/brian.a.mckinley/index/Sbicolor_313_v3.1.gene_exons.gtf -e -o ${filename}.temp.StringTie.gtf;
rm ${filename}.sam ${filename}.bam ${filename}.sorted.bam ${filename}.sh;
EOF
bsub < ${filename}.sh
done

for filename in *.temp.StringTie.gtf
do
awk '$3 == "transcript"' ${filename} > ${filename}.txt;
awk -F '[,;"]' '{print $17,$5,$16}' ${filename}.txt > ${filename}1.txt;
sort -k 2 ${filename}1.txt > ${filename}2.txt;
awk -F '[ ]' '{print $1}' ${filename}2.txt > ${filename}3.txt;
sed -i -e "1i ${filename}" ${filename}3.txt;
awk -F '[ ]' '{print $2}' ${filename}2.txt > transcriptIDs.txt;
done
awk -F '[ ]' '{print $1}' transcriptIDs.txt > transcriptIDs1.txt;
sed -i -e "1i TranscriptID" transcriptIDs1.txt;
paste -d ' ' *gtf3.txt > TPM.txt;
paste -d ' ' transcriptIDs1.txt TPM.txt > transcriptsTPM.txt;
rm *1.txt *2.txt *3.txt transcriptIDs.txt transcriptIDs1.txt TPM.txt;







#######################
## RNA-seq
#######################
for filename in /scratch/user/brian.a.mckinley/internode/*.fastq.gz
do
cat > ${filename}.sh << EOF
#BSUB -J starch3
#BSUB -L /bin/bash
#BSUB -W 04:00
#BSUB -n 5
#BSUB -R "span[ptile=5]"
#BSUB -R "select[nxt]"
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000
#BSUB -o starch2_stdout.%J
#BSUB -e starch2_stderr.%J
module load HISAT2/2.0.5-intel-2015B-Python-2.7.10
module load SAMtools/1.3-intel-2015B
module load StringTie/1.3.3b-intel-2015B
hisat2 -p 20 -q -x /scratch/user/brian.a.mckinley/index/sbi_base -U ${filename} -S ${filename}.sam --known-splicesite-infile /scratch/user/brian.a.mckinley/referanceV3/Sbicolorsplicesites.txt;
samtools view -@ 20 -bS ${filename}.sam > ${filename}.bam;
samtools sort -@ 20 ${filename}.bam -o ${filename}.sorted.bam;
stringtie ${filename}.sorted.bam -p 20 -G /scratch/user/brian.a.mckinley/index/Sbicolor_313_v3.1.gene_exons.gtf -e -o ${filename}.temp.StringTie.gtf;
rm ${filename}.sam ${filename}.bam ${filename}.sorted.bam ${filename}.sh;
EOF
bsub < ${filename}.sh
done


rm ${filename}.sam ${filename}.bam ${filename}.sorted.bam ${filename}.sh;


#stringtie ${filename}.sorted.bam -p 20 -e -o /scratch/user/brian.a.mckinley/test/testCTAB/Sample-"${filename##*/}"/Sample-"${filename##*/}.StringTie.gtf" -b /scratch/user/brian.a.mckinley/test/testCTAB/Sample-"${filename##*/}"/


################################
###StringTie Merge (this may not be useful because the code above specifies the use of a specific set of transcripts)
################################
#module load StringTie/1.3.1-intel-2015B
#ls -1 /scratch/user/qwerpoiu/Internode/*.temp.StringTie.gtf > /scratch/user/qwerpoiu/Internode/mergelist.txt
#stringtie --merge -p 8 -G /scratch/user/qwerpoiu/index/Sbicolor_313_v3.1.gene_exons.gtf -o /scratch/user/qwerpoiu/Internode/Internode_merged.gtf /scratch/user/qwerpoiu/Internode/mergelist.txt
################################
### combine all expression data from transcripts
################################
for filename in *.temp.StringTie.gtf
do
awk '$3 == "transcript"' ${filename} > ${filename}.txt
awk -F '[,;"]' '{print $17,$5,$16}' ${filename}.txt > ${filename}1.txt
sort -k 2 ${filename}1.txt > ${filename}2.txt
awk -F '[ ]' '{print $1}' ${filename}2.txt > ${filename}3.txt
sed -i -e "1i ${filename}" ${filename}3.txt
awk -F '[ ]' '{print $2}' ${filename}2.txt > transcriptIDs.txt
done
awk -F '[ ]' '{print $1}' transcriptIDs.txt > transcriptIDs1.txt
sed -i -e "1i TranscriptID" transcriptIDs1.txt
paste -d ' ' *gtf3.txt > TPM.txt
paste -d ' ' transcriptIDs1.txt TPM.txt > transcriptsTPM.txt
################################
### Erased all unnecessary files
################################
rm *.bam *.sam *.sh core* *Out *.StringTie.gtf *Abundance Internode_std* *csv *txt
###################################
###Run prepDE.py to extract counts-is there a way to extract TPM also?
###################################
module load python/2.7.10-intel-2015B
python /scratch/user/brian.a.mckinley/scripts/prepDE3.py -i /scratch/user/brian.a.mckinley/Della/DE/ -l 151 -t /scratch/user/brian.a.mckinley/Della/DE/counts
tr ',' ' ' < counts > counts.txt
(head -n 1 counts.txt && tail -n +3 counts.txt | sort) > counts_sorted.txt
paste -d ' ' transcriptsTPM.txt counts_sorted.txt > transcriptsTPM_counts.txt
##############################
###'R' script
##############################
#module load R_tamu/3.3.1-iomkl-2015B-default-mt
#R
# Load the following librairies, I am not sure if all of these are useful for our analysis
#library(ballgown)
#bg = ballgown(dataDir="/scratch/user/qwerpoiu/Internode/InternodeCTAB",samplePattern='Sample-', meas='all')
#whole_tx_table = *expr(bg, 'all')
#write.table(whole_tx_table, file ="InternodeAllData.csv",sep=",", col.names=TRUE)
##################################################
#   edgeR - calculation of Differential expression
##################################################
#biocLite("edgeR","limma","splines")
library(edgeR)
library(limma)
library(locfit)
library(splines)
# import data
#> setwd("/scratch/user/microchemmy")
raw.data <- read.csv("/scratch/user/qwerpoiu/edgeRInput.csv", header = TRUE, sep = ",")
# Specify which columns are associated with which sample groups
group = factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
# Create DGE list
y <- DGEList(counts=raw.data[, 1:18], genes=raw.data[, 0], group=group)
levels(y$samples$group)
# Filter gene set: keep only genes expresed greater than 0.5 counts per million in greater than 2 samples
keep <- rowSums(cpm(y) > 0.5) > 2
summary(keep)
y <- y[keep, ,keep.lib.sizes=FALSE]
# Perform TMM normalization
y <- calcNormFactors(y)
y$samples
#Perform multidimnesional scaling to visualize variation between biological replicates
plotMDS(y, top=500, main = "MDS of Count Data", cex=.75, labels = colnames(y$counts))
pdf ("MDS_plot_test.pdf" , width = 7 , height = 7) #inches
    plotMDS( y , main = "MDS Plot for Count Data", labels = colnames( y$counts) )
dev.off()
# Show TMM normalization factors
y$samples
# Plot histogram of library depth
barplot(y$samples$lib.size*1e-6, names=1:18, xpd = TRUE, axes = TRUE, col = 491, border = 300, xlab = "Samples", ylab="Reads aligned to exonic sequence (millions)")
# Create design matrix
design <- model.matrix(~0+group)
# insert GC bias estimation here refer to page 13 of edgeR manual
# Estimate common and tagwise dispersions in one function:
y <- estimateDisp(y,design)
# Mean-Variance plot to visualize how well the dispersion parameters fit the data
#meanVarPlot <- plotMeanVar( y , show.raw.vars=TRUE ,
#show.tagwise.vars=TRUE ,
#show.binned.common.disp.vars=FALSE ,
#show.ave.raw.vars=FALSE ,
#dispersion.method = "qcml" , NBline = TRUE ,
#nbins = 100 ,
#pch = 16 ,
#xlab ="Mean Expression (Log10 Scale)" ,
#ylab = "Variance (Log10 Scale)" ,
#main = "Mean-Variance Plot" )
# plot the biological coefficient of variation
#plotBCV(y)
# Qausi-likelihood dispersion estimation
fit <- glmQLFit(y,design)
#plotQLDisp(fit, cex = 0.2)
# Perform qausi-likelihood (QL) F-test to test for differentail expression
qlf <- glmQLFTest(fit, contrast = c(1,0,-1,0,0,0))
#To perform likelihood ratio tests:
#fit <- glmFit(y,design)
#lrt <- glmLRT(fit,coef=2)
#topTags(lrt)
# Print the top tags
lessFDR = sum(qlf$table$FDR <0.05)
result <- topTags(qlf, n="lessFDR")
write.table(result, file = "edgeROutput.txt")
options(max.print=1000000)
sink("test1.txt")
result
sink()



###########################
### Merge Wha tis the utility of this
############################
module add StringTie
ls -1 /scratch/user/qwerpoiu/Internode/*.temp.StringTie.gtf > /scratch/user/qwerpoiu/Internode/InternodeGTF_mergelist.txt
stringtie --merge -p 8 -G /scratch/user/qwerpoiu/index/Sbicolor_313_v3.1.gene_exons.gtf -o /scratch/user/qwerpoiu/Internode/Internode_stringTie_merged.gtf /scratch/user/qwerpoiu/Internode/InternodeGTF_mergelist.txt
