#!/bin/bash

NUMBER=0
SLEEP=300
FILE=output.txt

while ((NUMBER<600)); do
    clear
    head -n 20 $FILE
    echo ::::::::::::::::::::::::::::::::::
    tail -n 50 $FILE
    free -h
    ls -l $FILE
    sleep $SLEEP
    $NUMBER++
done
