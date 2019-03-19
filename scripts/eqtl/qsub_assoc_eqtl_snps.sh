#!/bin/bash
#
#PBS –N assocEqtlSnps
#PBS –M vpoten@gmail.com
#PBS -m abe
#PBS -l nodes=1:ppn=10,mem=8gb

java -Xmx8g -jar $PBS_O_HOME/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar assocEqtlSnps --props=$PBS_O_HOME/HlaQtl/ngsengine_eqtltrans.properties --output=$PBS_O_HOME/eqtl_trait_assoc --input=$PBS_O_HOME/eqtl_trans --eqtlDir=$PBS_O_HOME/eqtl_out --snps=$PBS_O_HOME/best_eqtls_1e05.raw.out --value=$PBS_O_HOME/GWAS_snps_diseases.csv > $PBS_O_HOME/assocEqtlSnps.out
