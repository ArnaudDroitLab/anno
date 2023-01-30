file=$1
echo $file
IFS=$'\n'
set -f
prefix=$(echo "$1" | awk -F '.raw_ref' '{print $1}')
echo "id,ensembl_gene,symbol,transcript_type" >> "$prefix.temp_no_alt_chr.csv"
if $(echo $file | grep -q "Ensembl"); then
        case $(echo $file | awk -F '.' ' { print $1 } ') in
        *"Bos_taurus")
                reg=":[1-2][0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Macaca_mulatta")
                reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Rattus_norvegicus")
                reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Mus_musculus")
                reg=":1[0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Homo_sapiens")
                reg=":1[0-9]:\|:2[0-2]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        esac
	write=true
	for i in $(cat $file); do
		if [[ $i == ">"* ]]; then
			if (( $(echo $i | grep -c $reg) == 0 )); then
				write=false
			else
				write=true
				echo $i | awk -F '.' ' { print $1 } ' >> "$prefix.no_alt_chr.fa"
				id=$(echo $i | awk -F '.' ' { print $1 } ' | cut -c 2-)
                        	ens=$(echo $i | awk -F 'gene:' ' { print $2 } ' | awk -F '.' '{print $1}')
                        	if (($(echo $i | grep -c 'gene_symbol') == 0)); then
                        	        sym="NA"
                        	else
                        	        sym=$(echo $i | awk -F 'gene_symbol:' ' { print $2 } ' | awk -F ' ' '{print $1}')
                        	fi
                        	tra=$(echo $i | awk -F 'transcript_biotype:' ' { print $2 } ' | awk -F ' ' '{print $1}')
				echo "$id,$ens,$sym,$tra" >> "$prefix.temp_no_alt_chr.csv"
			fi
		else
			if $write; then
				echo $i >> "$prefix.no_alt_chr.fa"
			fi
		fi
	done
else
	reg="PAR_Y"
	write=true
	for i in $(cat $file); do
                if [[ $i == ">"* ]]; then
                        if (( $(echo $i | grep -c $reg) == 1 )); then
                                write=false
                        else
                                write=true
                                echo $i | awk -F '|' ' { print $1 } ' | sed -e 's/\.[0-9]*//g' >> "$prefix.no_alt_chr.fa"
                		id=$(echo $i | awk -F '|' ' { print $1 } ' | sed -e 's/\.[0-9]*//g' | cut -c 2-)
                        	ens=$(echo $i | awk -F '|' ' { print $2 } ' | awk -F '.' '{print $1}')
                        	if [ $(echo $i | awk -F '|' ' { print $6 } ') == '-' ]; then
                        	        sym="NA"
                        	else
                        	        sym=$(echo $i | awk -F '|' ' { print $6 } ')
                        	fi
                        	tra=$(echo $i | awk -F '|' ' { print $8 } ')
				echo "$id,$ens,$sym,$tra" >> "$prefix.temp_no_alt_chr.csv"
			fi
		else
                        if $write; then
                                echo $i >> "$prefix.no_alt_chr.fa"
			fi
		fi
	done
fi
Rscript entrezid.R $(echo $prefix | awk -F '/' '{print $2}' | awk -F '.' '{print $1}') "$prefix.temp_no_alt_chr.csv" "$prefix.no_alt_chr.csv"
