filename1="A250.txt"
filename2="B250.txt"

cluster1="cluster.0_370.txt"
cluster2="cluster.1_93.txt"

rm -rf A0.txt
rm -rf A1.txt
while read -r line
do
    echo "$line"
    cat $cluster1 | grep -w "$line" >> A0.txt
    cat $cluster2 | grep -w "$line" >> A1.txt
done < "$filename1"

rm -rf B0.txt
rm -rf B1.txt
while read -r line
do
    echo "$line"
    cat $cluster1 | grep -w "$line" >> B0.txt
    cat $cluster2 | grep -w "$line" >> B1.txt
done < "$filename2"

