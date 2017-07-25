
# Example of DPR for genetic prediction

# Input File Format
DPR requires three input files containing genotypes, phenotypes and relatedness matrix. Genotype and phenotype files can be in two formats, either both in the PLINK binary ped format or both in the BIMBAM format. Mixing genotype and phenotype files from the two formats (for example, using PLINK files for genotypes and using BIMBAM files for phenotypes) will result in unwanted errors. BIMBAM format is particularly useful for imputed genotypes, as PLINK codes genotypes using 0/1/2, while BIMBAM can accommodate any real values between 0 and 2 (and any real values if paired with “-notsnp” option).

If you downloaded the DPR source code recently, you will find an “example” folder containing a small GWAS example dataset. This data set comes from the heterogeneous stock mice data, kindly provided by Wellcome Trust Centre for Human Genetics on the public domain http://mus.well.ox.ac.uk/mouse/HS/, with detailed described in [11]. The data set consists of 1,904 individuals from 85 families, all descended from eight inbred progenitor strains. We selected two phenotypes from this data set: the percentage of CD8+ cells, with measurements in 1,410 individuals; mean corpuscular hemoglobin (MCH), with measurements in 1,580 individuals. A total of 1,197 individuals have both phenotypes. The phenotypes were already corrected for sex, age, body weight, season and year effects by the original study, and we further quantile normalized the phenotypes to a standard normal distribution. In addition, we obtained a total of 12,226 autosomal SNPs, with missing genotypes replaced by the mean genotype of that SNP in the family. 

Genotype and phenotype files are in both BIMBAM and PLINK binary formats. For demonstration purpose, for CD8, we randomly divided the 85 families into two sets, where each set contained roughly half of the individuals (i.e. inter-family split) as in [6]. Therefore, the phenotype files contain four columns of phenotypes. The first column contains the quantitative phenotypes CD8 for all individuals. The second column contains quantitative phenotypes CD8 for individuals in the training set. The third column contains quantitative phenotypes CD8 for individuals in the test set. The fourth column contains the quantitative phenotypes MCH for all individuals.

# To fit a quantitative trait (i.e. CD8) using DPR with VB algorithm
./bin/DPR -g ./example/mouse_hs1940.geno.txt.gz -p ./example/mouse_hs1940.pheno.txt -n 2 -a ./example/mouse_hs1940.anno.txt -k  ./example/mouse_hs1940.cXX.txt -dpr 1 -nk 4 -o mouse_hs1940_CD8_vb

Explain

“-g” specifies BIMBAM genotypes, “-p” specifies phenotypes, “-a” specifies annotation file, “-k” specifies relatedness matrix, “-dpr 1” specifies fitting DPR using VB algorithm, “-nk 4” specifies four normal components included in into the mixture prior, “-o” specifies the output file.

./bin/DPR -bfile ./example/mouse_hs1940 -n 2 -k ./example/mouse_hs1940.cXX.txt -dpr 1 -nk 4 -o mouse_hs1940_CD8_vb

Explain

“-bfile” specifies plink files, i.e. mouse_hs1940.fam, mouse_hs1940.bim and mouse_hs1940.ped, “-n” specifies phenotypes using the 7th column of mouse_hs1940.fam, “-k” specifies relatedness matrix, “-dpr 1” specifies fitting DPR using VB algorithm, “-nk 4” specifies four normal components included in into the mixture prior, “-o” specifies the output file.


