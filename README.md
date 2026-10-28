# NetworkEvolution 
**Zhen-Zhu Xiao** 2026-10-28 17:15:02

### Downloaded data
Phenotype and genotype data for *Brassica rapa*: https://www.nature.com/articles/s41467-024-46841-4#Sec24, https://www.nature.com/articles/s41467-024-49383-x#data-availability, (NCBI sequence reads: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1105729) and https://onlinelibrary.wiley.com/doi/full/10.1111/mec.17811. 


### Network analysis

#### For each treatment network, I generate 1,000 random networks by rewiring edges while keeping each node's degree unchanged
code RandomRewire_M, Fig2_S4_TreatmentNetwork, 


#### Louvain algorithm and optimal algorith comparison 
I generate 10,000 networks with identical input to observe the distribution of modularity (*M*) with the Louvain algorithm. I found 100 runs of the Louvain algorithm are sufficient for the highest *M*. 

code: FigS3_ModuDistr

#### network simulation
Based on the observed edge changes between G1 and G10, I generate a randomized network by rewiring edges while preserving the same number of edges as in the original network.

code: Fig3

#### NMI between G1 and G10; NMI between treatment networks
code: Table1_FigS5_EdgeSimu_M_NMI_TrtDiff


### GWAS and PGS estimate

#### Bioinformatics
After variant call, I did hard filteration: 
    
    vcftools --gzvcf filtered_snps_onlypass.vcf.gz \
    --remove-indels \
    --min-alleles 2 --max-alleles 2 \
    --minGQ 30 \
    --minDP 5 --maxDP 100 \
    --min-meanDP 5 --max-meanDP 100 \
    --maf 0.01 \
    --max-missing 0.9 \
    --recode --recode-INFO-all \
    --out snps_filtered


#### GEMMA 

    # plink convert file format
    vcftools --gzvcf snps_renamed.vcf --plink --out all_snps
    plink --file all_snps --make-bed --out all_snps

    # generate binary ped files

    plink --file all_snps --make-bed --out all_snps

    gemma -bfile all_snps -gk 1 -o RelatedMatrix_center # generate relateness matrix 
    gemma -bfile all_snps -k RelatedMatrix_center.cXX.txt -n 6 -c ./covariate.cov -lmm 1 -o Trait11_cov

    code: FigS6_GWAS_Plot.R

- ##### Modular pleiotropy plot 
    code Fig4_Pleiotropy.R 

- ##### Concordant changes between phenotypic and genetic correlations: 
    code: FigS2_PearsonCorr.R


#### Polygenic score
    code: PGS_submit.sh



