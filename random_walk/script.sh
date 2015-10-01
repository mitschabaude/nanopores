#! /bin/bash
file=~/Desktop/script_test/file.txt
i=0
while [ $i -le 100 ]
do
    var="$(python python.py)"
    echo $var >> $file
    let i=$i+1
done
python plot_stat.py
