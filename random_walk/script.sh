#! /bin/bash
i=0
while [ $i -le 100 ]
do
    python run.py
    let i=$i+1
done
