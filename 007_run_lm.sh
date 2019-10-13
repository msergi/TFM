#!/usr/bin/env bash

for i in {1..22}
do
echo $i
Rscript 006_linear_m_subsets.R $i
done 2> x.log