#!/bin/bash

../a.out < input.txt > output.txt 2>error.txt
../a.out < input2.txt > output2.txt 2>error2.txt
for (( j=0; j<token; j++ ))
do
	mkdir $j
	cd $j
	echo "0.1 1 10.0 ${j} token" > input.txt
	cp ../../train.py .
	python3 train.py < input.txt > output.txt &

cd ..
done

wait
