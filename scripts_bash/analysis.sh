#class_begin=0
#class_end=22
#n_begin=1
#n_end=4

class_begin=0
class_end=0
n_begin=1
n_end=20
node=8

for ((i=$class_begin; i<=$class_end; i++)); do
    for ((j=$n_begin; j<=$n_end; j++)); do
        if [ "$i" = "0" ]; then
            DIR=modeling${j}
        else
            DIR=modeling${i}${j}
        fi

        # max ~25 secs per frame ~3456 frames per day ~48384 frames for two weeks
        echo ""; echo $DIR

        for ((k=0; k<$node; k++)); do
            process_output.py -f $DIR/pre-EM2D_output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_skDSS ISDCrossLinkMS_Data_Score_scEDC rmf_file rmf_frame_index
            #process_output.py -f $DIR/pre-EM2D_output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_scEDC rmf_file rmf_frame_index
            #process_output.py -f $DIR/pre-EM2D_output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_skDSS ISDCrossLinkMS_Data_Score_scEDC ISDCrossLinkMS_Data_Score_cliques rmf_file rmf_frame_index
        done

        # max ~100 secs per frame ~864 frames per day ~12096 frames for two weeks
        echo ""; echo $DIR

        for ((k=0; k<$node; k++)); do
            process_output.py -f $DIR/output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_skDSS ISDCrossLinkMS_Data_Score_scEDC rmf_file rmf_frame_index ElectronMicroscopy2D_None 
            #process_output.py -f $DIR/output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_scEDC rmf_file rmf_frame_index ElectronMicroscopy2D_None 
            #process_output.py -f $DIR/output/stat.$k.out -s Stopwatch_None_delta_seconds SimplifiedModel_Total_Score_None ExcludedVolumeSphere_None SimplifiedModel_Linker_Score_None ISDCrossLinkMS_Data_Score_scDSS ISDCrossLinkMS_Data_Score_skDSS ISDCrossLinkMS_Data_Score_scEDC ISDCrossLinkMS_Data_Score_cliques rmf_file rmf_frame_index ElectronMicroscopy2D_None 
        done

        sleep 1
    done
done

#MonteCarlo_Nframe
