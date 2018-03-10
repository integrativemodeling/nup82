SJ_EXTRACT_CHAINS_PY="./analysis_extract_chains2.py"
SJ_CHAIN_PY="./analysis_chain2.py"
SJ_COMBINE_CHAINS_PY="./analysis_combine_chains2.py"
#MODELLER="python"
MODELLER="mod9.15"

for FILES in *.pdb
do
    NEW_FILES=${FILES/.pdb/}

    ## Extract chains A-D separately
    python $SJ_EXTRACT_CHAINS_PY -pdb $FILES

    pulchra chainA.pdb; mv chainA.rebuilt.pdb ${NEW_FILES}_A.rebuilt.pdb; rm -rf chainA.pdb
    pulchra chainB.pdb; mv chainB.rebuilt.pdb ${NEW_FILES}_B.rebuilt.pdb; rm -rf chainB.pdb
    pulchra chainC.pdb; mv chainC.rebuilt.pdb ${NEW_FILES}_C.rebuilt.pdb; rm -rf chainC.pdb
    pulchra chainD.pdb; mv chainD.rebuilt.pdb ${NEW_FILES}_D.rebuilt.pdb; rm -rf chainD.pdb
    pulchra chainE.pdb; mv chainE.rebuilt.pdb ${NEW_FILES}_E.rebuilt.pdb; rm -rf chainE.pdb
    pulchra chainF.pdb; mv chainF.rebuilt.pdb ${NEW_FILES}_F.rebuilt.pdb; rm -rf chainF.pdb
    pulchra chainG.pdb; mv chainG.rebuilt.pdb ${NEW_FILES}_G.rebuilt.pdb; rm -rf chainG.pdb
    pulchra chainH.pdb; mv chainH.rebuilt.pdb ${NEW_FILES}_H.rebuilt.pdb; rm -rf chainH.pdb
    pulchra chainI.pdb; mv chainI.rebuilt.pdb ${NEW_FILES}_I.rebuilt.pdb; rm -rf chainI.pdb
    pulchra chainJ.pdb; mv chainJ.rebuilt.pdb ${NEW_FILES}_J.rebuilt.pdb; rm -rf chainJ.pdb
    pulchra chainK.pdb; mv chainK.rebuilt.pdb ${NEW_FILES}_K.rebuilt.pdb; rm -rf chainK.pdb
    pulchra chainL.pdb; mv chainL.rebuilt.pdb ${NEW_FILES}_L.rebuilt.pdb; rm -rf chainL.pdb
    pulchra chainM.pdb; mv chainM.rebuilt.pdb ${NEW_FILES}_M.rebuilt.pdb; rm -rf chainM.pdb
    pulchra chainN.pdb; mv chainN.rebuilt.pdb ${NEW_FILES}_N.rebuilt.pdb; rm -rf chainN.pdb
    pulchra chainO.pdb; mv chainO.rebuilt.pdb ${NEW_FILES}_O.rebuilt.pdb; rm -rf chainO.pdb
    pulchra chainP.pdb; mv chainP.rebuilt.pdb ${NEW_FILES}_P.rebuilt.pdb; rm -rf chainP.pdb
    pulchra chainQ.pdb; mv chainQ.rebuilt.pdb ${NEW_FILES}_Q.rebuilt.pdb; rm -rf chainQ.pdb

    : '
    ## Sort Ca atoms in the correct order of sequence
    grep -vwE "(END|TER)" chainA.pdb > chainA2.pdb
    sort -nk 6 chainA2.pdb > ${NEW_FILES}_A.pdb
    rm -rf chainA.pdb && rm -rf chainA2.pdb

    grep -vwE "(END|TER)" chainB.pdb > chainB2.pdb
    sort -nk 6 chainB2.pdb > ${NEW_FILES}_B.pdb
    rm -rf chainB.pdb && rm -rf chainB2.pdb

    grep -vwE "(END|TER)" chainC.pdb > chainC2.pdb
    sort -nk 6 chainC2.pdb > ${NEW_FILES}_C.pdb
    rm -rf chainC.pdb && rm -rf chainC2.pdb

    grep -vwE "(END|TER)" chainD.pdb > chainD2.pdb
    sort -nk 6 chainD2.pdb > ${NEW_FILES}_D.pdb
    rm -rf chainD.pdb && rm -rf chainD2.pdb

    grep -vwE "(END|TER)" chainE.pdb > chainE2.pdb
    sort -nk 6 chainE2.pdb > ${NEW_FILES}_E.pdb
    rm -rf chainE.pdb && rm -rf chainE2.pdb

    grep -vwE "(END|TER)" chainF.pdb > chainF2.pdb
    sort -nk 6 chainF2.pdb > ${NEW_FILES}_F.pdb
    rm -rf chainF.pdb && rm -rf chainF2.pdb

    grep -vwE "(END|TER)" chainG.pdb > chainG2.pdb
    sort -nk 6 chainG2.pdb > ${NEW_FILES}_G.pdb
    rm -rf chainG.pdb && rm -rf chainG2.pdb

    grep -vwE "(END|TER)" chainH.pdb > chainH2.pdb
    sort -nk 6 chainH2.pdb > ${NEW_FILES}_H.pdb
    rm -rf chainH.pdb && rm -rf chainH2.pdb

    grep -vwE "(END|TER)" chainI.pdb > chainI2.pdb
    sort -nk 6 chainI2.pdb > ${NEW_FILES}_I.pdb
    rm -rf chainI.pdb && rm -rf chainI2.pdb

    grep -vwE "(END|TER)" chainJ.pdb > chainJ2.pdb
    sort -nk 6 chainJ2.pdb > ${NEW_FILES}_J.pdb
    rm -rf chainJ.pdb && rm -rf chainJ2.pdb

    #mv chainC.pdb ${NEW_FILES}_C.pdb
    #mv chainD.pdb ${NEW_FILES}_D.pdb
   
    ## Reconstruct side chain atoms
    pulchra ${NEW_FILES}_A.pdb
    pulchra ${NEW_FILES}_B.pdb
    pulchra ${NEW_FILES}_C.pdb
    pulchra ${NEW_FILES}_D.pdb
    pulchra ${NEW_FILES}_E.pdb
    pulchra ${NEW_FILES}_F.pdb
    pulchra ${NEW_FILES}_G.pdb
    pulchra ${NEW_FILES}_H.pdb
    pulchra ${NEW_FILES}_I.pdb
    pulchra ${NEW_FILES}_J.pdb
    
    rm -rf ${NEW_FILES}_A.pdb && rm -rf ${NEW_FILES}_B.pdb
    rm -rf ${NEW_FILES}_C.pdb && rm -rf ${NEW_FILES}_D.pdb
    rm -rf ${NEW_FILES}_E.pdb && rm -rf ${NEW_FILES}_F.pdb
    rm -rf ${NEW_FILES}_G.pdb && rm -rf ${NEW_FILES}_H.pdb
    rm -rf ${NEW_FILES}_I.pdb && rm -rf ${NEW_FILES}_J.pdb
    '
done

# recover chain ID
$MODELLER $SJ_CHAIN_PY

# combine chains A-D into a single PDB file
for FILES in *_A.rebuilt.pdb
do
    NEW_FILES=${FILES/_A.rebuilt.pdb/}

    python $SJ_COMBINE_CHAINS_PY -pdb1 $FILES -pdb2 ${NEW_FILES}_B.rebuilt.pdb -pdb3 ${NEW_FILES}_C.rebuilt.pdb -pdb4 ${NEW_FILES}_D.rebuilt.pdb -pdb5 ${NEW_FILES}_E.rebuilt.pdb -pdb6 ${NEW_FILES}_F.rebuilt.pdb -pdb7 ${NEW_FILES}_G.rebuilt.pdb -pdb8 ${NEW_FILES}_H.rebuilt.pdb -pdb9 ${NEW_FILES}_I.rebuilt.pdb -pdb10 ${NEW_FILES}_J.rebuilt.pdb -pdb11 ${NEW_FILES}_K.rebuilt.pdb -pdb12 ${NEW_FILES}_L.rebuilt.pdb -pdb13 ${NEW_FILES}_M.rebuilt.pdb -pdb14 ${NEW_FILES}_N.rebuilt.pdb -pdb15 ${NEW_FILES}_O.rebuilt.pdb -pdb16 ${NEW_FILES}_P.rebuilt.pdb -pdb17 ${NEW_FILES}_Q.rebuilt.pdb -out ${NEW_FILES}a.pdb
    grep -vwE "(END|ENDMDL|MODEL      0)" ${NEW_FILES}a.pdb > ${NEW_FILES}r.pdb

    rm -rf $FILES && rm -rf ${NEW_FILES}_B.rebuilt.pdb
    rm -rf ${NEW_FILES}_C.rebuilt.pdb && rm -rf ${NEW_FILES}_D.rebuilt.pdb
    rm -rf ${NEW_FILES}_E.rebuilt.pdb && rm -rf ${NEW_FILES}_F.rebuilt.pdb
    rm -rf ${NEW_FILES}_G.rebuilt.pdb && rm -rf ${NEW_FILES}_H.rebuilt.pdb
    rm -rf ${NEW_FILES}_I.rebuilt.pdb && rm -rf ${NEW_FILES}_J.rebuilt.pdb
    rm -rf ${NEW_FILES}_K.rebuilt.pdb && rm -rf ${NEW_FILES}_L.rebuilt.pdb
    rm -rf ${NEW_FILES}_M.rebuilt.pdb && rm -rf ${NEW_FILES}_N.rebuilt.pdb
    rm -rf ${NEW_FILES}_O.rebuilt.pdb && rm -rf ${NEW_FILES}_P.rebuilt.pdb
    rm -rf ${NEW_FILES}_Q.rebuilt.pdb
    rm -rf ${NEW_FILES}a.pdb

done
