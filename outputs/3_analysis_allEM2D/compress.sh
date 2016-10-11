PWD_PARENT=$(pwd)

#for DIR_SA in modeling_*
#do
#    cd $DIR_SA

    #for DIR in modeling161_failed modeling174_failed; do
    #for DIR in modeling*_failed; do
    for DIR in kmeans_*; do
        if [ -d "$DIR" ]; then
            echo "tar czf $DIR.tar.gz $DIR"
            tar czf $DIR.tar.gz $DIR
            rm -rf $DIR
            #mv $DIR.tar.gz back_up
        fi
    done
    exit -1
    
    for DIR in analysis*; do
        if [ -d "$DIR" ]; then
            echo "tar czf $DIR.tar.gz $DIR"
            tar czf $DIR.tar.gz $DIR
            rm -rf $DIR
        fi
    done


#    cd ..
#done

-rw-r--r--. 1 sjkim sali 0 Sep 23 08:32 n82_c10.po9133593.104
-rw-r--r--. 1 sjkim sali 0 Sep 23 19:23 n82_c11b.po9133596.116
-rw-r--r--. 1 sjkim sali 0 Sep 23 06:58 n82_c11b.po9133596.118
-rw-r--r--. 1 sjkim sali 0 Sep 23 14:54 n82_c13b.po9133600.136
-rw-r--r--. 1 sjkim sali 0 Sep 23 15:52 n82_c13b.po9133600.139
-rw-r--r--. 1 sjkim sali 0 Sep 23 07:42 n82_c14.po9133601.143
-rw-r--r--. 1 sjkim sali 0 Sep 23 13:58 n82_c16.po9133605.161
-rw-r--r--. 1 sjkim sali 0 Sep 23 12:36 n82_c5.po9133583.51
-rw-r--r--. 1 sjkim sali 0 Sep 23 13:15 n82_c5.po9133583.53
-rw-r--r--. 1 sjkim sali 0 Sep 23 15:30 n82_c6b.po9133586.67
-rw-r--r--. 1 sjkim sali 0 Sep 23 10:35 n82_c6b.po9133586.68
-rw-r--r--. 1 sjkim sali 0 Sep 23 21:37 n82_c7b.po9133588.78
-rw-r--r--. 1 sjkim sali 0 Sep 23 16:08 n82_c7b.po9133588.79
-rw-r--r--. 1 sjkim sali 0 Sep 23 09:37 n82_c7.po9133587.72
-rw-r--r--. 1 sjkim sali 0 Sep 23 15:33 n82_c8.po9133589.83
