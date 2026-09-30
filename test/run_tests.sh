#!/bin/bash

# Smoke tests for clipplotr using the bundled CD55 test data
# Runs in a temporary directory so nothing is written to the repository
# Usage: test/run_tests.sh [output directory to keep plots]

set -u

TEST_DIR=$(cd "$(dirname "$0")" && pwd)
CLIPPLOTR="$TEST_DIR/../clipplotr"
WORK_DIR=$(mktemp -d)
OUT_DIR=${1:-$WORK_DIR/plots}
mkdir -p "$OUT_DIR"
OUT_DIR=$(cd "$OUT_DIR" && pwd)
trap 'rm -rf "$WORK_DIR"' EXIT

cp "$TEST_DIR"/*.gz "$TEST_DIR"/*.bigwig "$TEST_DIR"/*.tsv "$WORK_DIR"/
cd "$WORK_DIR"

XLINKS='test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz,test_hnRNPC_iCLIP_rep2_LUjh25_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_ctrl_rep1_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_ctrl_rep2_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_KD1_rep2_all_xlink_events.bedgraph.gz,test_U2AF65_iCLIP_KD2_rep1_all_xlink_events.bedgraph.gz'
LABELS='hnRNPC 1,hnRNPC 2,U2AF2 WT 1,U2AF2 WT 2,U2AF2 KD 1,U2AF2 KD 2'
CLIP_GROUPS='hnRNPC,hnRNPC,U2AF2 WT,U2AF2 WT,U2AF2 KD,U2AF2 KD'
SIZE_FACTORS='4.869687,9.488133,1.781117,10.135903,4.384385,8.227587'
COVERAGE='test_ERR127306_plus.bigwig,test_ERR127307_plus.bigwig,test_ERR127302_plus.bigwig,test_ERR127303_plus.bigwig'
COVERAGE_GROUPS='CTRL,CTRL,KD,KD'
GTF='CD55_gencode.v34lift37.annotation.gtf.gz'
REGION='chr1:207513000:207515000:+'

PASSED=0
FAILED=0

# expect_pass <name> <args...>: must exit 0 and write <name>.png
expect_pass() {
  local name=$1; shift
  if Rscript "$CLIPPLOTR" -g "$GTF" --cache_dir "$WORK_DIR/cache" -o "$OUT_DIR/$name.png" "$@" > "$WORK_DIR/$name.log" 2>&1 && [ -s "$OUT_DIR/$name.png" ]; then
    echo "PASS  $name"; PASSED=$((PASSED + 1))
  else
    echo "FAIL  $name"; FAILED=$((FAILED + 1)); tail -5 "$WORK_DIR/$name.log" | sed 's/^/      /'
  fi
}

# expect_fail <name> <expected message> <args...>: must exit non-zero with the message
expect_fail() {
  local name=$1 expected=$2; shift 2
  if Rscript "$CLIPPLOTR" -g "$GTF" --cache_dir "$WORK_DIR/cache" -o "$OUT_DIR/$name.png" "$@" > "$WORK_DIR/$name.log" 2>&1; then
    echo "FAIL  $name (exited 0)"; FAILED=$((FAILED + 1))
  elif grep -q -- "$expected" "$WORK_DIR/$name.log"; then
    echo "PASS  $name"; PASSED=$((PASSED + 1))
  else
    echo "FAIL  $name (message not found: $expected)"; FAILED=$((FAILED + 1)); tail -5 "$WORK_DIR/$name.log" | sed 's/^/      /'
  fi
}

# Full example from the README
expect_pass readme -x "$XLINKS" -l "$LABELS" --groups "$CLIP_GROUPS" -n custom --size_factors "$SIZE_FACTORS" \
  -s rollmean -w 50 -y test_Alu_rev.bed.gz --auxiliary_labels 'reverse Alu' \
  --coverage "$COVERAGE" --coverage_groups "$COVERAGE_GROUPS" -r "$REGION" --highlight 207513650:207513800 -a transcript

# Option combinations
for norm in none libsize maxpeak libsize_maxpeak custom custom_maxpeak; do
  for smooth in none rollmean; do
    expect_pass "norm_${norm}_${smooth}" -x "$XLINKS" -l "$LABELS" --groups "$CLIP_GROUPS" --size_factors "$SIZE_FACTORS" \
      -n "$norm" -s "$smooth" -w 50 -r "$REGION" -a none
  done
done

for annot in transcript gene none; do
  expect_pass "annot_${annot}" -x "$XLINKS" -l "$LABELS" -r "$REGION" -a "$annot" -y test_Alu_rev.bed.gz --coverage "$COVERAGE"
  expect_pass "annot_${annot}_flip" -x "$XLINKS" -l "$LABELS" -r "$REGION" -a "$annot" -y test_Alu_rev.bed.gz --coverage "$COVERAGE" --flip_x
done

# Region given as gene name and gene id (with and without version)
expect_pass region_gene_name -x "$XLINKS" -r CD55 -a transcript
expect_pass region_gene_id -x "$XLINKS" -r ENSG00000196352.16_8 -a gene
expect_pass region_gene_id_unversioned -x "$XLINKS" -r ENSG00000196352 -a gene

# Default labels are unique even when file names share a prefix
expect_pass default_labels -x "$XLINKS" -r "$REGION" -a none

# Capitalised .bedGraph extension
cp test_hnRNPC_iCLIP_rep1_LUjh03_all_xlink_events.bedgraph.gz hnRNPC.bedGraph.gz
expect_pass bedGraph_extension -x hnRNPC.bedGraph.gz -r "$REGION" -a none

# Auxiliary BED with no features in the region
expect_pass auxiliary_empty -x "$XLINKS" -r chr1:207514500:207515000:+ -w 20 -y test_Alu_rev.bed.gz -a transcript

# Small region
expect_pass small_region -x "$XLINKS" -r chr1:207513700:207513750:+ -w 10 -a transcript

# Opposite strand: no transcripts in the region
expect_pass no_transcripts -x "$XLINKS" -r chr1:207513000:207515000:- -a transcript

# Errors exit non-zero with a helpful message
expect_fail unknown_gene "was not found" -x "$XLINKS" -r NOTAGENE
expect_fail bad_region "should be given as" -x "$XLINKS" -r chr1:100:200
expect_fail labels_length "has 2 entries" -x "$XLINKS" -l 'a,b' -r "$REGION"
expect_fail duplicate_labels "must be unique" -x "$XLINKS" -l 'a,a,b,c,d,e' -r "$REGION"
expect_fail size_factors_length "size_factors has" -x "$XLINKS" -n custom --size_factors 1,2 -r "$REGION"
expect_fail gaussian_removed "gaussian smoothing has been removed" -x "$XLINKS" -s gaussian -r "$REGION"
expect_fail original_removed "original annotation style has been removed" -x "$XLINKS" -a original -r "$REGION"
expect_fail missing_file "do not exist" -x missing.bedgraph -r "$REGION"
expect_fail window_too_large "larger than the region" -x "$XLINKS" -r chr1:207513700:207513750:+ -w 100

# Samplesheet instead of comma-separated options (same plot as the README example)
expect_pass samplesheet --samples test_samples.tsv -n custom -s rollmean -w 50 -r "$REGION" --highlight 207513650:207513800 -a transcript
mkdir -p sheet_dir && sed 's/\ttest_/\t..\/test_/' test_samples.tsv > sheet_dir/samples.tsv
expect_pass samplesheet_relative_paths --samples sheet_dir/samples.tsv -n custom -r "$REGION" -a none
expect_fail samplesheet_with_xlinks "cannot be combined" --samples test_samples.tsv -x "$XLINKS" -r "$REGION"

# Several regions: {region} in the output name gives one file each
printf 'CD55\n# a comment\nchr1:207513500:207514000:+\n' > regions.txt
if Rscript "$CLIPPLOTR" -g "$GTF" --cache_dir "$WORK_DIR/cache" -x "$XLINKS" -r chr1:207513000:207515000:+ --regions regions.txt \
     -w 50 -o "$OUT_DIR/multi_{region}.png" > multi.log 2>&1 \
   && [ -s "$OUT_DIR/multi_chr1_207513000_207515000_+.png" ] && [ -s "$OUT_DIR/multi_CD55.png" ] && [ -s "$OUT_DIR/multi_chr1_207513500_207514000_+.png" ]; then
  echo "PASS  multi_region_files"; PASSED=$((PASSED + 1))
else
  echo "FAIL  multi_region_files"; FAILED=$((FAILED + 1)); tail -5 multi.log | sed 's/^/      /'
fi

# Several regions in one .pdf: one page each
if Rscript "$CLIPPLOTR" -g "$GTF" --cache_dir "$WORK_DIR/cache" -x "$XLINKS" -r 'CD55,chr1:207513500:207514000:+' -w 50 -o "$OUT_DIR/multi.pdf" > multi_pdf.log 2>&1 \
   && [ "$(grep -ac '/Type /Page\b' "$OUT_DIR/multi.pdf")" -eq 2 ]; then
  echo "PASS  multi_region_pdf"; PASSED=$((PASSED + 1))
else
  echo "FAIL  multi_region_pdf"; FAILED=$((FAILED + 1)); tail -5 multi_pdf.log | sed 's/^/      /'
fi

expect_fail multi_region_png "needs to contain {region}" -x "$XLINKS" -r 'CD55,chr1:207513500:207514000:+'

echo
echo "$PASSED passed, $FAILED failed"
[ "$FAILED" -eq 0 ]
