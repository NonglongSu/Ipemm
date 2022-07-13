echo "bash runbash100.sh 100 500"
imax=$1
jmax=$2
for ((i=1;i<=${imax};i++))
do
  subD=Gs/$i
  echo "Working on Subfolder-$i"
  for ((j=1;j<=${jmax};j++))
  do
    bash gap_rm.sh ${subD}/${j}.fa Gs_trim/$i/${j}.fasta
  done
done
