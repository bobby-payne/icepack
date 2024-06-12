#!/bin/bash

for file in *.tar; do
	echo "Extracting $file"
	tar -xf "$file"
	rm $file
done
