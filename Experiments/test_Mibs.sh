#!/bin/bash
#test script

data_folder="/home/xus6/BFKO/Dataset/DNegMibs/"
result_file="results.log"
echo $data_folder

file="DNeg"

rm mibs.log
rm $result_file

for ((n=10; n<=50; n+=5)) do
  for ((r=0;r<10;r++)) do
    for ((j=0;j<10;j++))do
      #file_name="${file}_n$n"
      file_name="${file}_n${n}k${r}j${j}"
      echo $file_name
      
      echo "*******************************************************************************************" >> mibs.log
      echo $file_name >> mibs.log
      echo "*******************************************************************************************" >> mibs.log
      
      start_time="$(date -u +%s.%N)"
      mibs -Alps_instance "$data_folder${file_name}.mps" -MibS_auxiliaryInfoFile "$data_folder${file_name}.txt" -Alps_timeLimit 600 >> mibs.log
      end_time="$(date -u +%s.%N)"
      elapsed="$(bc <<<"$end_time-$start_time")"
      
      printf "$file_name \t $elapsed \n" >> $result_file
    done
  
  done
done
