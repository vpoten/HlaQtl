#!/bin/bash
#
#PBS –N eqtlTransFilter
#PBS –M vpoten@gmail.com
#PBS -m abe
#PBS -l mem=8gb

java -Xmx8g -jar $PBS_O_HOME/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar eqtlTransFilter --props=$PBS_O_HOME/HlaQtl/ngsengine_eqtltrans.properties --output=$PBS_O_HOME/best_eqtl_trans_stats --input=$PBS_O_HOME/eqtl_trans --blast=$PBS_O_HOME/software/ncbi-blast-2.2.28+/bin --eqtlDir=$PBS_O_HOME/eqtl_out --snps=$PBS_O_HOME/best_eqtls_1e05.raw.out > $PBS_O_HOME/eqtlTransFilter.out


