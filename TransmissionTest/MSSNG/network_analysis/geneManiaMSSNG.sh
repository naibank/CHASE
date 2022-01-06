#!/bin/bash -x

#PBS -N gm.mssng
#PBS -l nodes=1:ppn=32
#PBS -l vmem=24g
#PBS -d /hpf/largeprojects/tcagstor/users/worrawat/CHASE/TranmissionTest/ShaniaOneDrive/MSSNG/network_analysis/

module load java

java -jar /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/GeneMania.jar QueryRunner --data /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/gmdata-2014-08-12-core/ MSSNG_sanitized.txt --out scores
