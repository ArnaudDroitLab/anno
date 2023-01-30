file=$1
if [[ $file =~ \.gz$ ]]; then
  	gzip -dk $file
 	file=${file:0:-3}
fi
echo "$file"
# basic minimal file (With 5 protein coding sequences)
if [[ "$file" == *"Gencode"* ]]; then
	elements=$(grep ">" $file | awk -F '|' ' { print $8 } ' | sort | uniq)
else
	elements=$(grep ">" $file | awk -F ' ' ' { print $5 } ' | sort | uniq)
fi

for line in $elements; do
	if [[ "$line" == *"gene_biotype:protein_coding"* || "$line" == "protein_coding" ]]; then
		line_n=$(grep -n "|$line|\| $line " $file | head -n 1 | awk -F ':' ' { print $1 } ')
        	line_n2=$(($(tail -n +$line_n $file | grep -n ">" | tail -n +2 | head -n 1 | awk -F ':' ' { print $1 } ')-1))
        	tail -n +"$line_n" $file | head -n $line_n2 >> "out/$file"
		for i in 1 2 3 4
		do
			prev_line=$line_n
			line_n=$(($(tail -n +$(($prev_line+$line_n2)) $file | grep -n "|$line|\| $line " | head -n 1 | awk -F ':' ' { print $1} ')+$prev_line+$line_n2-1))
        		line_n2=$(($(tail -n +$line_n $file | grep -n ">" | tail -n +2 | head -n 1 | awk -F ':' ' { print $1 } ')-1))
			tail -n +$line_n $file | head -n $line_n2 >> "out/$file"
		done
	else
		line_n=$(grep -n "|$line|\| $line " $file | head -n 1 | awk -F ':' ' { print $1 } ')
		line_n2=$(($(tail -n +$line_n $file | grep -n ">" | tail -n +2 | head -n 1 | awk -F ':' ' { print $1 } ')-1))
		tail -n +"$line_n" $file | head -n $line_n2 >> "out/$file"
	fi
done

# Add a non standard sequence for no alt chr step
if $(echo $file | grep -q "Ensembl"); then
	case $(echo $file | awk -F '.' ' { print $1 } ') in
	"Bos_taurus")
		reg=":[1-2][0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
	"Macaca_mulatta")
		reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
	"Rattus_norvegicus")
		reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
	"Mus_musculus")
		reg=":1[0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
	"Homo_sapiens")
		reg=":1[0-9]:\|:2[0-2]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
	esac

	if (( $(grep ">" "out/$file" | grep -v -c $reg ) == 0 && $(grep ">" $file | grep -v -c $reg) != 0 )); then
		line_n=$(grep -n ">" $file | grep -v $reg | head -n 1 | awk -F ':' ' { print $1 } ')
		line_n2=$(($(tail -n +$line_n $file | grep -n ">" | tail -n +2 | head -n 1 | awk -F ':' ' { print $1 } ')-1))
		tail -n +"$line_n" $file | head -n $line_n2 >> "out/$file"
	fi
elif $(echo $file | grep -q "Homo_sapiens.Gencode"); then
	reg="PAR_Y"
	if ! (( $(grep -c -q $reg "out/$file") )); then
		line_n=$(grep -n $reg $file | head -n 1 | awk -F ':' ' { print $1 } ')
                line_n2=$(($(tail -n +$line_n $file | grep -n ">" | tail -n +2 | head -n 1 | awk -F ':' ' { print $1 } ')-1))
		tail -n +"$line_n" $file | head -n $line_n2 >> "out/$file"
	fi
fi
