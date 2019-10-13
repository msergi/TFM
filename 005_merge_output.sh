#!/usr/bin/env bash

for i in {1..22}
do
o="$o output$i.txt "
done
cat $o > all_output.txt