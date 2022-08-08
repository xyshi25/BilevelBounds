#!/bin/bash
# The script for instances solved by Mibs

##########################################################################
#
## run experiments for the bilevel knapsack interdiction problem
#
###########################################################################

data_folder="Dataset/DNegMibs/"
result_file="DNeg_results_mibs.log"
log_path="Log/Log_KIP"
echo $data_folder

file="DNeg"

rm -f $result_file

mkdir -p Log/Log_KIP



for ((n=10; n<=50; n+=10)) do
  for ((r=0;r<10;r++)) do
    for ((j=0;j<10;j++))do
      #file_name="${file}_n$n"
      file_name="${file}_n${n}k${r}j${j}"
      echo "$data_folder${file_name}"
      
      
      start_time="$(date -u +%s.%N)"
      mibs -Alps_instance "$data_folder${file_name}.mps" -MibS_auxiliaryInfoFile "$data_folder${file_name}.txt" -Alps_timeLimit 600 > "${log_path}/${file_name}_mibs.log"
      end_time="$(date -u +%s.%N)"
      elapsed="$(bc <<<"$end_time-$start_time")"
      
      printf "$file_name \t $elapsed \n" >> $result_file
    done
  
  done
done




##########################################################################
#
## run experiments for the bilevel vertex cover problem
#
###########################################################################

## sysmmetric objective

data_folder="Dataset/BVCMibs/"
result_file="BVC_Symm_results_mibs.log"
log_path="Log/Log_BVC_Symm"
echo $data_folder

file="BVC"

rm -f $result_file

mkdir -p Log/Log_BVC_Symm


for ((n=10; n<=50; n+=5)) do
  deg=$((($n+1)/2))
  for ((j=0;j<10;j++))do
    file_name="${file}_n${n}d${deg}b${deg}j${j}"
    
    echo "$data_folder${file_name}"
    start_time="$(date -u +%s.%N)"
    mibs -Alps_instance "$data_folder${file_name}.mps" -MibS_auxiliaryInfoFile "$data_folder${file_name}.txt" -Alps_timeLimit 1800 > "${log_path}/${file_name}_mibs.log"
    end_time="$(date -u +%s.%N)"
    elapsed="$(bc <<<"$end_time-$start_time")"
    
    printf "$file_name \t $elapsed \n" >> $result_file
    
    
  done
  
done

# asymmetric objective 
data_folder="Dataset/BVC_ASYMM_Mibs/"
result_file="BVC_ASymm_results_mibs.log"
log_path="Log/Log_BVC_ASymm"
echo $data_folder

file="BVC"

rm -f $result_file

mkdir -p Log/Log_BVC_ASymm


for ((n=10; n<=50; n+=5)) do
  deg=$((($n+1)/2))
  for ((j=0;j<10;j++))do
    file_name="${file}_n${n}d${deg}b${deg}j${j}"
    
    echo "$data_folder${file_name}"
    start_time="$(date -u +%s.%N)"
    mibs -Alps_instance "$data_folder${file_name}.mps" -MibS_auxiliaryInfoFile "$data_folder${file_name}.txt" -Alps_timeLimit 1800 > "${log_path}/${file_name}_mibs.log"
    end_time="$(date -u +%s.%N)"
    elapsed="$(bc <<<"$end_time-$start_time")"
    
    printf "$file_name \t $elapsed \n" >> $result_file
    
    
  done
  
done
