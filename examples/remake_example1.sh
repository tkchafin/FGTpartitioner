#!/bin/bash

rm example_1.vcf.gz*
cp example_noHet.vcf example_1.vcf
bgzip example_1.vcf
tabix -h -f -p vcf example_1.vcf.gz
