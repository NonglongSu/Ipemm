infile=$1
oufile=$2
sed  's/-//g' ${infile} > ${oufile}  
