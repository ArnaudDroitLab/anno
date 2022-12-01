file=$1
echo $file
IFS=$'\n'
set -f
prefix=$(echo "$1" | awk -F '.raw_ref' '{print $1}')
echo "id,ensembl_gene,symbol,transcript_type" >> "$prefix.temp_protein_coding.csv"
for i in $(cat $file); do
	if [[ $i == ">"* ]]; then
                if [[ "$i" == *" gene_biotype:protein_coding "* || "$i" == *"|protein_coding|"* ]]; then
                	write=true
			if [[ $file == *"Gencode"* ]]; then
                        	echo $i | awk -F '|' ' { print $1 } ' | sed -e 's/\.[0-9]*//g' >> "$prefix.protein_coding.fa"
				id=$(echo $i | awk -F '|' ' { print $1 } ' | sed -e 's/\.[0-9]*//g' | cut -c 2-)
                        	ens=$(echo $i | awk -F '|' ' { print $2 } ' | awk -F '.' '{print $1}')
                        	if [ $(echo $i | awk -F '|' ' { print $6 } ') == '-' ]; then
                        	        sym="NA"
                        	else
                        	        sym=$(echo $i | awk -F '|' ' { print $6 } ')
                        	fi
                        	tra=$(echo $i | awk -F '|' ' { print $8 } ')
			else
                		echo $i | awk -F '.' ' { print $1 } ' >> "$prefix.protein_coding.fa"
				id=$(echo $i | awk -F '.' ' { print $1 } ' | cut -c 2-)
                        	ens=$(echo $i | awk -F 'gene:' ' { print $2 } ' | awk -F '.' '{print $1}')
                        	if (($(echo $i | grep -c 'gene_symbol') == 0)); then
                        	        sym="NA"
                        	else
                        	        sym=$(echo $i | awk -F 'gene_symbol:' ' { print $2 } ' | awk -F ' ' '{print $1}')
                        	fi
                        	tra=$(echo $i | awk -F 'transcript_biotype:' ' { print $2 } ' | awk -F ' ' '{print $1}')
                	fi
			echo "$id,$ens,$sym,$tra" >> "$prefix.temp_protein_coding.csv"
		else
			write=false
		fi
        else
                if $write; then
                       echo $i >> "$prefix.protein_coding.fa"
                fi
        fi
done

Rscript entrezid.R $(echo $prefix | awk -F '/' '{print $2}' | awk -F '.' '{print $1}') "$prefix.temp_protein_coding.csv" "$prefix.protein_coding.csv"
