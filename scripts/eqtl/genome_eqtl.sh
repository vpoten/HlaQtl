#!/bin/bash

jars=/home/users/ipb/HlaQtl/HlaQtl-1.0-SNAPSHOT-bin/dist/lib/NGSUtils-1.0-SNAPSHOT.jar
script=/home/users/ipb/HlaQtl/genome_eqtl.groovy

groovy -cp $jars $script $1

