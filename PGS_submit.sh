#!/bin/bash 
#SBATCH --job-name=presice28
#SBATCH --mail-user=zhenzhux@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10G
#SBATCH --time=3:30:00
#SBATCH --account=sigbio_project19
#SBATCH --partition=sigbio
#SBATCH --output=/home/zhenzhux/log/%x-%j.out
#SBATCH --error=/home/zhenzhux/log/%x-%j.err

conda init bash
conda activate new
cd /scratch/sigbio_project_root/sigbio_project19/zhenzhux/gwas/GPScore/polyScore42

######## all trait all dir loop 
for dir in /scratch/sigbio_project_root/sigbio_project19/zhenzhux/gwas/GPScore/*_gwas
do
    echo "Processing directory: ${dir}"

    cd "${dir}" || continue

    prefix=$(basename "${dir}" _gwas)

    mkdir -p LD0.2_result

    for BASE in Trait*_sumStat.txt.gz
    do
        TRAIT=$(basename "${BASE}" _sumStat.txt.gz)

        echo "Running PRSice for ${TRAIT} in ${prefix}"

        Rscript /scratch/sigbio_project_root/sigbio_project19/zhenzhux/gwas/PRSice/PRSice.R \
            --prsice /scratch/sigbio_project_root/sigbio_project19/zhenzhux/gwas/GPScore/PRSice/PRSice_linux \
            --base "${BASE}" \
            --target "${prefix}_target" \
            --binary-target F \
            --pheno "${prefix}_target_traits.pheno" \
            --pheno-col "${TRAIT}" \
            --score std \
            --num-auto 10 \
            --cov "${prefix}_target_cov.eigenvec" \
            --keep-ambig \
            --clump-kb 30 \
            --clump-r2 0.2 \
            --base-maf MAF:0.01 \
            --a1 A1 \
            --a2 A2 \
            --stat Beta \
            --thread 5 \
            --out "./LD0.2_result/${TRAIT}_${prefix}_LD0.2_PGS"

    done
done








