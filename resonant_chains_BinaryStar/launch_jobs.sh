#!/bin/bash

for i in {0..75}
do
   start=$((i * 2))
   end=$(((i+1) * 2))
   echo "Slice $star:$end"
   nohup python resonance_chain_planetary_system.py $start $end >text/out${i}.txt &
done
