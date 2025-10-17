samp_id=${2}
cell_threshold_count=${3}

for i in `ls $1`; do
    echo ${i}
    tail -n +2 ${i}  |\
    cut -d '	' -f 2-  |\
    sed 's/FALSE/0/g' |\
    sed 's/TRUE/1/g' |\
    sed 's/	//g' |\
    sort | uniq -c |\
    sort -n |\
    sed 's/\s\+/ /g' |\
    grep -v '"' |\
    cut -f 2- |\
    awk -v ctc=${cell_threshold_count} '$1 > ctc' |\
    sed 's/ /	/g' |\
    cut -f 2- > intermediates/$samp_id.vertices_size.tsv
done

for i in `ls $1`; do
    echo ${i}
    tail -n +2 ${i}  |\
    cut -d '	' -f 2-   |\
    sed 's/FALSE/0/g' |\
    sed 's/TRUE/1/g' |\
    sed 's/	//g' |\
    sort | uniq -c |\
    sort -n |\
    sed 's/\s\+/ /g' |\
    grep -v '"' |\
    cut -f 2- |\
    awk -v ctc=${cell_threshold_count} '$1 > ctc' |\
    sed 's/ /	/g' |\
    cut -f 3 > intermediates/$samp_id.vertices.tsv
done

python scripts/hamming.py -vi intermediates/$samp_id.vertices.tsv -si intermediates/$samp_id.vertices_size.tsv -el intermediates/$samp_id.edgelist -so intermediates/$samp_id.edgelist_meta -ci $1 -vco intermediates/$samp_id.node_pmhcs -cco intermediates/$samp_id.pmhc_counts
