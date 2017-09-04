
# Example of DPR for genetic prediction

## Input File Format
DPR requires three input files containing genotypes, phenotypes and relatedness matrix. Genotype and phenotype files can be in two formats, either both in the PLINK binary ped format or both in the BIMBAM format. Mixing genotype and phenotype files from the two formats (for example, using PLINK files for genotypes and using BIMBAM files for phenotypes) will result in unwanted errors. BIMBAM format is particularly useful for imputed genotypes, as PLINK codes genotypes using 0/1/2, while BIMBAM can accommodate any real values between 0 and 2 (and any real values if paired with “-notsnp” option).

If you downloaded the DPR source code recently, you will find an “example” folder containing a small GWAS example dataset. This data set comes from the heterogeneous stock mice data, kindly provided by Wellcome Trust Centre for Human Genetics on the public domain [here](http://mus.well.ox.ac.uk/mouse/HS/). The data set consists of 1,904 individuals from 85 families, all descended from eight inbred progenitor strains. We selected two phenotypes from this data set: the percentage of CD8+ cells, with measurements in 1,410 individuals; mean corpuscular hemoglobin (MCH), with measurements in 1,580 individuals. A total of 1,197 individuals have both phenotypes. The phenotypes were already corrected for sex, age, body weight, season and year effects by the original study, and we further quantile normalized the phenotypes to a standard normal distribution. In addition, we obtained a total of 12,226 autosomal SNPs, with missing genotypes replaced by the mean genotype of that SNP in the family. 

Genotype and phenotype files are in both BIMBAM and PLINK binary formats. For demonstration purpose, for CD8, we randomly divided the 85 families into two sets, where each set contained roughly half of the individuals (i.e. inter-family split). Therefore, the phenotype files contain four columns of phenotypes. The first column contains the quantitative phenotypes CD8 for all individuals. The second column contains quantitative phenotypes CD8 for individuals in the training set. The third column contains quantitative phenotypes CD8 for individuals in the test set. The fourth column contains the quantitative phenotypes MCH for all individuals.

## Some examples
### To fit a quantitative trait (i.e. CD8) using DPR with VB algorithm
***`./bin/DPR -g ./example/mouse_hs1940.geno.txt.gz -p ./example/mouse_hs1940.pheno.txt -n 2 -a ./example/mouse_hs1940.anno.txt -k  ./example/mouse_hs1940.cXX.txt -dpr 1 -nk 4 -o mouse_hs1940_CD8_vb`***

**Explain:**

“-g” specifies BIMBAM genotypes, “-p” specifies phenotypes, “-a” specifies annotation file, “-k” specifies relatedness matrix, “-dpr 1” specifies fitting DPR using VB algorithm, “-nk 4” specifies four normal components included in into the mixture prior, “-o” specifies the output file.

***`./bin/DPR -bfile ./example/mouse_hs1940 -n 2 -k ./example/mouse_hs1940.cXX.txt -dpr 1 -nk 4 -o mouse_hs1940_CD8_vb`***

**Explain:**

“-bfile” specifies plink files, i.e. mouse_hs1940.fam, mouse_hs1940.bim and mouse_hs1940.ped, “-n” specifies phenotypes using the 7th column of mouse_hs1940.fam, “-k” specifies relatedness matrix, “-dpr 1” specifies fitting DPR using VB algorithm, “-nk 4” specifies four normal components included in into the mixture prior, “-o” specifies the output file.

### To fit a quantitative trait (i.e. CD8) using DPR with MCMC algorithm
***`./bin/DPR -bfile ./example/mouse_hs1940 -n 2 -k ./example/mouse_hs1940.cXX.txt -dpr 2 -nk 4 -m 100000 -t 1 -w 10000 -s 10000 -o mouse_hs1940_CD8_mcmc`***

**Explain:**

“-bfile” specifies plink files, i.e. mouse_hs1940.fam, mouse_hs1940.bim and mouse_hs1940.ped, “-n” specifies phenotypes using the 7th column of mouse_hs1940.fam, “-k” specifies relatedness matrix, “-dpr 2” specifies fitting DPR using MCMC algorithm, “-nk 4” specifies four normal components included in into the mixture prior, "-m 100000" specifies the top 100000 SNPs to be included in the model (if the value is larger than the total number of the SNPs, then all variants are used), "-t 1" update those non-selected, likely unimportant SNPs once every 1 iteration, “-w 10000” specifies 10000 burn-ins, “-s 10000” specifies 10000 samplings after burn-in, “-o” specifies the output file.

### To fit a quantitative trait (i.e. CD8) using adaptive LDR
***`./bin/DPR -bfile ./example/mouse_hs1940 -n 2 -k ./example/mouse_hs1940.cXX.txt -dpr 3 -mnk 6 -sp 0.2 -w 10000 -s 10000 -o mouse_hs1940_CD8_ada`***

**Explain:**

“-bfile” specifies plink files, i.e. mouse_hs1940.fam, mouse_hs1940.bim and mouse_hs1940.ped, “-n” specifies phenotypes using the 7th column of mouse_hs1940.fam, “-k” specifies relatedness matrix, “-dpr 3” specifies fitting DPR using adaptive LDR, “-mnk 6” specifies the maximum number of normal components included in the mixture prior, “-sp 0.2” specifies samplings for this adaptive process, “-w 10000” specifies 10000 burn-ins, “-s 10000” specifies 10000 samplings after burn-in, “-o” specifies the output file. More specifically, this command will repeat fitting MCMC sampling for nk=2 to nk=6, for each nk (i.e. 2, 3, 4, 5 and 6), the burn-in number is 2000=10000×0.2 (w×sp), the MC sampling number is also 2000=10000×0.2 (s×sp). After the optimal nk (say nk*) is selected based on the smallest DIC across nk=2 to nk=6, this command continuous perform DPR using MCMC sampling with nk* normal components, and now the burn-in number is 10000, the MC sampling number is 10000.

## Cite
#### Ping Zeng and Xiang Zhou (2017). Non-parametric genetic prediction of complex traits with latent Dirichlet process regression models. Nature Communications [in press](http://www.biorxiv.org/content/early/2017/06/13/149609).

## Reference
1. Ferguson TS. A Bayesian analysis of some nonparametric problems. Ann Stat. 1973; 1: 209-230.
2. Andrews DF, Mallows CL. Scale mixtures of normal distributions. J R Stat Soc Ser B. 1974; 36: 99-102.
3. Sethuraman J. A constructive definition of Dirichlet priors. Stat Sinica. 1994; 4: 639 - 650.
4. Ishwaran H, James LF. Gibbs sampling methods for stick-breaking priors. J Am Stat Assoc. 2001; 96.
5. Blei DM, Jordan MI. Variational inference for Dirichlet process mixtures. Bayesian Anal. 2006; 1: 121-143.
6. Zhou X, Carbonetto P, Stephens M. Polygenic modeling with Bayesian sparse linear mixed models. PLoS Genet. 2013; 9: e1003264.
7. Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies. Nat Genet. 2012; 44: 821-824.
8. Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MA, Bender D, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007; 81: 559-575.
9. Guan Y, Stephens M. Practical Issues in Imputation-Based Association Mapping. PLoS Genet. 2008; 4: e1000279.
10. Howie BN, Donnelly P, Marchini J. A Flexible and Accurate Genotype Imputation Method for the Next Generation of Genome-Wide Association Studies. PLoS Genet. 2009; 5: e1000529.
11. Valdar W, Solberg LC, Gauguier D, Burnett S, Klenerman P, Cookson WO, et al. Genome-wide genetic association of complex traits in heterogeneous stock mice. Nat Genet. 2006; 38: 879-887.


