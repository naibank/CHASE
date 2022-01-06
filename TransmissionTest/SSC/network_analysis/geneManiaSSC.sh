#!/bin/bash -x

#PBS -N gm.ssc
#PBS -l nodes=1:ppn=32
#PBS -l vmem=24g
#PBS -d /hpf/largeprojects/tcagstor/users/worrawat/CHASE/TranmissionTest/ShaniaOneDrive/SSC/network_analysis/

module load java

java -jar /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/GeneMania.jar QueryRunner --data /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/gmdata-2014-08-12-core/ SSC_probands_sanitized.txt --out scores
java -jar /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/GeneMania.jar QueryRunner --data /hpf/largeprojects/tcagstor/users/worrawat/outdate_projects/gene-mania/gmdata-2014-08-12-core/ SSC_unaffectedsiblings_sanitized.txt --out scores
