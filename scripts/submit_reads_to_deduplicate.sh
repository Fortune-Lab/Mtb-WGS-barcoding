
while IFS=$'\t' read -r -a array;
do
	if [ -f "${array[0]}.samtools.coverage" ]; then
		echo "${array[0]} completed."
	else	
		sbatch Remove_Duplicates.sh ${array[0]} ${array[1]} ${array[2]}
		sleep 1
	fi
done < $1
