
# Example of DPR for genetic prediction

# Input File Format
DPR requires three input files containing genotypes, phenotypes and relatedness matrix. Genotype and phenotype files can be in two formats, either both in the PLINK binary ped format or both in the BIMBAM format. Mixing genotype and phenotype files from the two formats (for example, using PLINK files for genotypes and using BIMBAM files for phenotypes) will result in unwanted errors. BIMBAM format is particularly useful for imputed genotypes, as PLINK codes genotypes using 0/1/2, while BIMBAM can accommodate any real values between 0 and 2 (and any real values if paired with “-notsnp” option).
