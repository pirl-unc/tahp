threshold=${1}
infile=${2}
outfile=${3}
sample_id=${4}

awk -v samp_id="$sample_id" -v thresh="$threshold" -v OFS="\t" '
NR==1 { ncol = NF }
{
    for (i=2; i<=NF; i++) {
        if ($i == "TRUE") count[i]++
        data[NR,i] = $i
    }
    data[NR,1] = $1
}
END {
    thresh_count = NR * thresh / 100
    print thresh_count > "intermediates/" samp_id ".cell_thresh_count"
    keep[1] = 1
    for (i=2; i<=ncol; i++) {
        if (count[i] >= thresh_count) keep[i] = 1
    }
    for (r=1; r<=NR; r++) {
        first = 1
        for (i=1; i<=ncol; i++) {
            if (keep[i]) {
                if (!first) printf OFS
                printf "%s", data[r,i]
                first = 0
            }
        }
        print ""
    }
}' "$infile" > "$outfile"
