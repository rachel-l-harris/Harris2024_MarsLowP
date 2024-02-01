#!/bin/bash

#SBATCH -J array_MarsGas_MAC_12mbar_blastn_ncbi-all-CDS
#SBATCH -n 1 #--ntasks
#SBATCH -c 16 #--cpus_per-task
#SBATCH -N 1 #node count
#SBATCH -t 0-4:00:00 #--run time in hours
#SBATCH -p shared #--partition
#SBATCH --mem=1G #--memory request in GB
#SBATCH -o /n/holyscratch01/girguis_lab/rachel_harris/Mars_low-pressure/MarsGas/MAC_12mbar/cds_KEGG-mby-all-CDS/mappedReads/blastn/array_MarsGas_MAC_12mbar_blastn_ncbi-all-CDS_mapping_%A_%a.out
#SBATCH -e /n/holyscratch01/girguis_lab/rachel_harris/Mars_low-pressure/MarsGas/MAC_12mbar/cds_KEGG-mby-all-CDS/mappedReads/blastn/array_MarsGas_MAC_12mbar_blastn_ncbi-all-CDS_mapping_%A_%a.err
#SBATCH --array=0-3
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=rachel_harris@fas.harvard.edu

echo "$SLURM_ARRAY_TASK_ID" #this will be 0,1,2,3

echo "cd /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS"
cd /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS

echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "####        makeblastdb       ###"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

echo ""
echo "Make a blast database of ncbi-all-CDS"
echo ""

#the PATH to makeblastdb is written into my .bash_profile so no need to activate a conda environment

#echo "makeblastdb -in ncbi.M.barkeri-MS_cds.fasta -dbtype nucl -parse_seqids -out blastdb/blastdb_ncbi.M.barkeri-MS_cds"
#makeblastdb -in ncbi.M.barkeri-MS_cds.fasta -dbtype nucl -parse_seqids -out blastdb/blastdb_ncbi.M.barkeri-MS_cds
#echo ""
#echo "blastdb generated"
#echo ""

echo "cd /n/holyscratch01/girguis_lab/rachel_harris/Mars_low-pressure"
cd /n/holyscratch01/girguis_lab/rachel_harris/Mars_low-pressure

echo ""
echo "####################################"
echo "BLASTn MarsGas/MAC_12mbar MAPPED READS"
echo "#####################################"
echo ""
echo "cd MarsGas/MAC_12mbar/cds_KEGG-mby-all-CDS/mappedReads/mappers"
cd MarsGas/MAC_12mbar/cds_KEGG-mby-all-CDS/mappedReads/mappers
echo ""

echo "Load the sample names into a bash array"
echo "samples=(12a_A1-MAC1 12a_A1-MAC2 12b_A1-MAC1 12b_A1-MAC2)"
samples=(12a_A1-MAC1 12a_A1-MAC2 12b_A1-MAC1 12b_A1-MAC2)
echo ""

echo "@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@"
echo "###    seqtk      ###"
echo "@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@"

#echo "source activate /n/holylfs05/LABS/girguis_lab/Lab/envs/seqtk/"
#source activate /n/holylfs05/LABS/girguis_lab/Lab/envs/seqtk/
#echo ""

#echo "convert fastq to fasta"
#echo "forward reads"
#echo "seqtk seq -a mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fastq > fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fasta"
#seqtk seq -a mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fastq > fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fasta
#echo "reverse reads"
#echo "seqtk seq -a mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fastq > fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fasta"
#seqtk seq -a mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fastq > fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fasta
#echo ""

echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "####          blastn           ###"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

#the PATH to blastn is written into my .bash_profile so no need to activate a conda environment

#echo "blastn -db /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS/blastdb/blastdb_ncbi.M.barkeri-MS_cds -query fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fasta -num_threads 16 -outfmt 11 -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.asn"
#blastn -db /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS/blastdb/blastdb_ncbi.M.barkeri-MS_cds -query fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.fasta -num_threads 16 -outfmt 11 -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.asn
#echo "blastn -db /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS/blastdb/blastdb_ncbi.M.barkeri-MS_cds -query fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fasta -num_threads 16 -outfmt 11 -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.asn"
#blastn -db /n/girguis_lab/Users/rachel_harris/Mars_low-pressure/M.barkeri_MS/blastdb/blastdb_ncbi.M.barkeri-MS_cds -query fasta_mappers/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.fasta -num_threads 16 -outfmt 11 -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.asn
#echo ""
#echo "blastn complete; formatting"
#echo "blast_formatter -archive ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.asn -outfmt "5 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out"
#blast_formatter -archive ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.asn -outfmt "5 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out
#echo "blast_formatter -archive ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.asn -outfmt "5 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out"
#blast_formatter -archive ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.asn -outfmt "5 qseqid qlen sseqid slen qstart qend sstart send length pident nident mismatch evalue bitscore staxids saccver stitle" -out ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out
#echo "complete"
#echo ""


echo ""
echo ""
echo ""
echo "BLASTn COMPLETE!"

echo "python /n/girguis_lab/Everyone/blast-QC/BLAST_QC_PYTHON/BLAST-QC.py -f ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out -ff XML -t n -o ../blastn/blast-qc/QC_mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out -n 1 -or e -e 0.001 -b 40 -i 0.9"
python /n/girguis_lab/Everyone/blast-QC/BLAST_QC_PYTHON/BLAST-QC.py -f ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out -ff XML -t n -o ../blastn/blast-qc/QC_mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R1_KEGG-mby-all-CDS.BLASTn.out -n 1 -or e -e 0.001 -b 40 -i 0.9

echo "python /n/girguis_lab/Everyone/blast-QC/BLAST_QC_PYTHON/BLAST-QC.py -f ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out -ff XML -t n -o ../blastn/blast-qc/QC_mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out -n 1 -or e -e 0.001 -b 40 -i 0.9"
python /n/girguis_lab/Everyone/blast-QC/BLAST_QC_PYTHON/BLAST-QC.py -f ../blastn/mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out -ff XML -t n -o ../blastn/blast-qc/QC_mappers_${samples[$SLURM_ARRAY_TASK_ID]}_R2_KEGG-mby-all-CDS.BLASTn.out -n 1 -or e -e 0.001 -b 40 -i 0.9

echo ""
echo "blast-qc complete"
