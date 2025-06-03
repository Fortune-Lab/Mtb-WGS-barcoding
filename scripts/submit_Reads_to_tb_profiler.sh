
while IFS=$'\t' read -r -a array;
do
	if [ -f "${array[0]}.tbprofiler.done" ]; then
                echo "${array[0]} completed."
        else
		sbatch run_tb-profiler.sh ${array[0]} ${array[1]} ${array[2]}
		sleep 1
	fi
done < $1
