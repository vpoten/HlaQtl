#!/bin/bash
#
#PBS –N controlGroupLD
#PBS –M vpoten@gmail.com
#PBS -m abe
#PBS -l nodes=1:ppn=10,mem=8gb

java -Xmx8g -jar $PBS_O_HOME/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar controlGroupLD --props=$PBS_O_HOME/HlaQtl/ngsengine_eqtltrans.properties --output=$PBS_O_HOME/gtex_immunochip_ld --eqtlDir=$PBS_O_HOME/eqtl_out --snps=$PBS_O_HOME/gtex_immunochip_ld/groups_gtex.all.txt > $PBS_O_HOME/gtex_immunochip_ld/controlGroupLD.out


