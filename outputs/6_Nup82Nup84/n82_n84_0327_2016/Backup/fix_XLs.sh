for FILE in *.pdb
do
    file=${FILE/.pdb/}
    rm -rf ${file}r.pdb
    echo $file

    while IFS= read -r line; do
        NEW_LINE=${line/CA  BEA C 624/CA  BEA C 622}
        NEW_LINE=${NEW_LINE/CA  BEA C 674/CA  BEA C 671}
        NEW_LINE=${NEW_LINE/CA  BEA D 624/CA  BEA D 622}
        NEW_LINE=${NEW_LINE/CA  BEA D 674/CA  BEA D 671}
        NEW_LINE=${NEW_LINE/CA  BEA L   6/CA  BEA L   2}
        NEW_LINE=${NEW_LINE/CA  BEA L 233/CA  BEA L 231}
        NEW_LINE=${NEW_LINE/CA  BEA J   6/CA  BEA J   2}
        NEW_LINE=${NEW_LINE/CA  BEA J 233/CA  BEA J 231}
        NEW_LINE=${NEW_LINE/CA  BEA O 167/CA  BEA O 168}
        NEW_LINE=${NEW_LINE/CA  BEA P 199/CA  BEA P 198}

        printf '%s\n' "$NEW_LINE" >> ${file}r.pdb
    done < "$file.pdb"
done
