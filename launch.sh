#!/bin/bash

export OMP_NUM_THREADS=18

for i in {1..2}
do
	mkdir ${i}
	cd ${i}

	echo "Debug 125000 250000000 25000 20000 10.0 train ${i} 0.05 0.013 200 Exit" > input.txt
	echo "Debug 250000000 380150000 25000 20000 10.0 test ${i} 0.05 0.013 200 Exit" > input2.txt

	cat ../job.sh | sed "s/token/${i}/g" > ./job.sh
	chmod +x ./job.sh
	./job.sh &

	cd ..
done

wait

grep "" */*/output.txt

