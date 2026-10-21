#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=siwinski@uvm.edu

###change to where your data is
cd /gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/population_genomics/mydata

###this is where you should output your data 
OUTDIR=/gpfs1/home/s/i/siwinski/projects/eco_genomics_2026/population_genomics/myresults

echo "This is my out directory"
echo $OUTDIR
###CallVCF

VCF=Lottia_gigantea_pl_flt3.renamed.vcf.gz

echo "This is my VCF"
echo $VCF
###load plink
module load plink

echo "This is my version of Plink"
echo plink --version

##LD thin first - try changing the windows and LD threshold?

plink --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.2 \
    --out $OUTDIR/Lg_LD_0.2_alldata

##now lets run the PCA - what other parameters can you add when making a pca in plink?

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract $OUTDIR/Lg_LD_0.2_alldata.prune.in \
    --make-bed --pca --out $OUTDIR/Lg_LD_0.2_alldata

echo "my script is done!"