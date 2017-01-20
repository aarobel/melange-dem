DT=0.002
TF=$((3600*24*5))
NT=`echo "$TF/$DT" | bc -l`

STR=500000
LFC=0.03

module load matlab

for HSI in 1; do

	STRL=$((${STR}/100000))
	LFCS=$(echo "scale=0; $LFC*1000/1" | bc) 
	HSIS=$(echo "scale=0; $HSI*10/1" | bc)
	matlab -nodesktop -nosplash -r "makebatch_nlf(${NT},0,'fjordbond_final_nlf',${HSI},${STR},${LFC});quit;"

	for ITER in $(seq 1 1 14); do
   		echo "${ITER}"
   		export ITER
   		LGHTSTR="sbatch -d singleton batch_fjordbond_final_nlf_h${HSIS}_str${STRL}e5_lfc${LFCS}.sh ${NT} ${ITER}"
   		eval $LGHTSTR
	done
done
