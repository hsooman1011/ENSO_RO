#!/bin/bash

### Starting and ending random_num to be newly created ###
start=1
end=9

### base RO code with random_num = 0 ###
base_RO="RO_sweep_000.py"



for (( i=$start; i<=$end; i++ ))
do
    ### Update the random_num parameter in the newly generated RO_sweep_{random_num}.py
    random_num=$(printf "%03d" $i)
    new_RO="RO_sweep_${random_num}.py"
	cp "$base_RO" "${new_RO}"
	if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    	sed -i '' "s/random_num = [0-9]*/random_num = $i/" "${new_RO}"
	else
    # Linux and others
    	sed -i "s/random_num = [0-9]*/random_num = $i/" "${new_RO}"
	fi

done

