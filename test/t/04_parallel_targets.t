#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
source "$root/test/test-simple-bash/lib/test-simple.bash" tests 2
caller=${FREEBAYES:-freebayes}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
ref="$root/test/tiny/q.fa"
bam="$root/test/tiny/NA12878.chr22.tiny.bam"
printf 'q\t0\t1000\nq\t1000\t2000\n' > "$work/targets.bed"
for n in 1 12; do
    "$caller" -f "$ref" -t "$work/targets.bed" --threads "$n" "$bam" > "$work/$n.vcf"
    sed '/^#/d' "$work/$n.vcf" > "$work/$n.records"
done
awk '$2 > 1000 {found=1} END {exit !found}' "$work/1.records"
ok cmp -s "$work/1.records" "$work/12.records" 'separate BED targets retain all serial records'
# Duplicate the reference and alignments under a second contig name.
cat "$ref" > "$work/ref.fa"
sed 's/^>q$/>q2/' "$ref" >> "$work/ref.fa"
samtools faidx "$work/ref.fa"
samtools view -H "$bam" | awk 'BEGIN {FS=OFS="\t"} {print; if ($1=="@SQ" && $2=="SN:q") { $2="SN:q2"; print }}' > "$work/two.sam"
samtools view "$bam" >> "$work/two.sam"
samtools view "$bam" | awk 'BEGIN {FS=OFS="\t"} {$1=$1"_q2"; $3="q2"; if ($7=="q") $7="q2"; print}' >> "$work/two.sam"
samtools view -b -o "$work/two.bam" "$work/two.sam"
samtools index "$work/two.bam"
for n in 1 12; do
    "$caller" -f "$work/ref.fa" --threads "$n" "$work/two.bam" > "$work/$n.vcf"
    sed '/^#/d' "$work/$n.vcf" > "$work/$n.records"
done
grep -q $'^q2\t' "$work/1.records"
ok cmp -s "$work/1.records" "$work/12.records" 'separate contigs retain all serial records'
