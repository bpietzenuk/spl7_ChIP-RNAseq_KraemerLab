#!/bin/bash
export LC_ALL=C

for file in $@
do
	sed '1d' $file > $file.noHeader.txt
done
