file=$1
echo $file
IFS=$'\n'
set -f
prefix=$(echo "$1" | awk -F '.raw_ref' '{print $1}')
if $(echo $file | grep -q "Ensembl"); then
        case $(echo $file | awk -F '.' ' { print $1 } ') in
        *"Bt")
                reg=":[1-2][0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Mmu")
                reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Rn")
                reg=":1[0-9]:\|:20:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Mm")
                reg=":1[0-9]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        *"Hs")
                reg=":1[0-9]:\|:2[0-2]:\|:[1-9]:\|:X:\|:Y:\|:MT:";;
        esac
        write=true
        for i in $(cat $file); do
                if [[ $i == ">"* ]]; then
                        if (( $(echo $i | grep -c $reg) == 0 )); then
                                write=false
                        else
                                write=true
                                echo $i >> "$prefix.no_alt_full.fa"
                        fi
                else
                        if $write; then
                                echo $i >> "$prefix.no_alt_full.fa"
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
                                echo $i >> "$prefix.no_alt_full.fa"
                        fi
                else
                        if $write; then
                                echo $i >> "$prefix.no_alt_full.fa"
                        fi
                fi
        done
fi
