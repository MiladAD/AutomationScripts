#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: pipeline.sh <vcfFile> <outputPrefix> <refFASTA> <bed> <annovarFolder> <annovardbFolder> <refBuild> <thread> [protocol] [operation]

  vcfFile:         Input VCF file (optionally bgzipped)
  outputPrefix:    Prefix for output files (can include output dir)
  refFASTA:        Reference FASTA file (indexed for bcftools norm)
  bed:             BED file (regions to KEEP with bcftools filter -R)
  annovarFolder:   Path to ANNOVAR scripts (contains table_annovar.pl)
  annovardbFolder: Path to ANNOVAR databases
  refBuild:        Reference build (e.g., hg19 or hg38)
  thread:          ANNOVAR threads
  protocol:        (Optional) ANNOVAR protocol string
  operation:       (Optional) ANNOVAR operation string

Notes:
  - This pipeline preserves samples/GT by running ANNOVAR on a sites-only VCF
    and then merging ANNOVAR INFO back into the original multi-sample VCF.
  - It also normalizes illegal INFO tag names emitted by some ANNOVAR/dbNSFP
    fields (e.g. "GERP++", "Eigen-", "M-CAP", "1000g...").
EOF
  exit 1
}

if [[ $# -lt 8 ]]; then
  echo "Error: Missing arguments." >&2
  usage
fi

vcfFile=$1
outputPrefix=$2
refFASTA=$3
bed=$4
annovarFolder=$5
annovardbFolder=$6
refBuild=$7
thread=$8

# Default protocol/operation (hg19 defaults; override with args 9/10)
protocol="refGeneWithVer,knownGene,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad211_genome,gnomad211_exome,exac03,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
operation="g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"

# hg38 example (use via args 9/10 or swap defaults above):
# protocol="refGeneWithVer,knownGene,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad41_exome,gnomad41_genome,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
# operation="g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"

if [[ $# -ge 10 ]]; then
  protocol=$9
  operation=${10}
fi

mkdir -p "$(dirname "${outputPrefix}")"

need_cmd() {
  local c="$1"
  command -v "$c" >/dev/null 2>&1 || { echo "Error: Required command not found: $c" >&2; exit 1; }
}

need_file() {
  local f="$1"
  [[ -e "$f" ]] || { echo "Error: File not found: $f" >&2; exit 1; }
}

need_cmd bcftools
need_cmd tabix
need_cmd bgzip
need_cmd sed
[[ -x "${annovarFolder}/table_annovar.pl" ]] || { echo "Error: ANNOVAR script not found/executable: ${annovarFolder}/table_annovar.pl" >&2; exit 1; }

need_file "$vcfFile"
need_file "$refFASTA"
need_file "$bed"
need_file "$annovardbFolder"

# Sed in-place portability (GNU vs BSD/macOS)
sed_inplace() {
  if sed --version >/dev/null 2>&1; then
    sed -i "$@"
  else
    sed -i '' "$@"
  fi
}

# Quick write permission check for output directory (catches "Permission denied" early)
assert_writable_dir() {
  local dir="$1"
  local tmp="${dir}/.write_test.$$"
  if ! touch "$tmp" 2>/dev/null; then
    echo "Error: cannot write to output directory: $dir" >&2
    echo "Fix permissions/ownership or choose a different outputPrefix path." >&2
    exit 1
  fi
  rm -f "$tmp"
}

# Fix invalid INFO tag names in both header (##INFO=<ID=...>) and body (INFO keys)
fix_vcf_info_tags_inplace() {
  local vcf="$1"

  # IMPORTANT:
  # - VCF tag IDs must not start with digits and must not contain + or -
  # - We patch known problematic tags from ANNOVAR/dbNSFP/1000g conventions.

  sed_inplace \
    \
    # 1000g2015aug_all starts with digit -> prefix with 'g'
    -e '/^##INFO=<ID=1000g2015aug_all([,>])/{s/ID=1000g2015aug_all/ID=g1000g2015aug_all/;}' \
    -e 's/;1000g2015aug_all=/;g1000g2015aug_all=/g' \
    \
    # "#Chr" invalid due to '#'
    -e '/^##INFO=<ID=#Chr([,>])/{s/ID=#Chr/ID=Chr/;}' \
    -e 's/;#Chr=/;Chr=/g' \
    \
    # Hyphenated tags -> underscores
    -e '/^##INFO=<ID=M-CAP_pred([,>])/{s/ID=M-CAP_pred/ID=M_CAP_pred/;}' \
    -e 's/;M-CAP_pred=/;M_CAP_pred=/g' \
    \
    -e '/^##INFO=<ID=LIST-S2_pred([,>])/{s/ID=LIST-S2_pred/ID=LIST_S2_pred/;}' \
    -e 's/;LIST-S2_pred=/;LIST_S2_pred=/g' \
    \
    -e '/^##INFO=<ID=fathmm-MKL_coding_pred([,>])/{s/ID=fathmm-MKL_coding_pred/ID=fathmm_MKL_coding_pred/;}' \
    -e 's/;fathmm-MKL_coding_pred=/;fathmm_MKL_coding_pred=/g' \
    \
    -e '/^##INFO=<ID=fathmm-XF_coding_pred([,>])/{s/ID=fathmm-XF_coding_pred/ID=fathmm_XF_coding_pred/;}' \
    -e 's/;fathmm-XF_coding_pred=/;fathmm_XF_coding_pred=/g' \
    \
    -e '/^##INFO=<ID=Eigen-raw_coding([,>])/{s/ID=Eigen-raw_coding/ID=Eigen_raw_coding/;}' \
    -e 's/;Eigen-raw_coding=/;Eigen_raw_coding=/g' \
    \
    -e '/^##INFO=<ID=Eigen-phred_coding([,>])/{s/ID=Eigen-phred_coding/ID=Eigen_phred_coding/;}' \
    -e 's/;Eigen-phred_coding=/;Eigen_phred_coding=/g' \
    \
    -e '/^##INFO=<ID=Eigen-PC-raw_coding([,>])/{s/ID=Eigen-PC-raw_coding/ID=Eigen_PC_raw_coding/;}' \
    -e 's/;Eigen-PC-raw_coding=/;Eigen_PC_raw_coding=/g' \
    \
    -e '/^##INFO=<ID=Eigen-PC-phred_coding([,>])/{s/ID=Eigen-PC-phred_coding/ID=Eigen_PC_phred_coding/;}' \
    -e 's/;Eigen-PC-phred_coding=/;Eigen_PC_phred_coding=/g' \
    \
    # GERP++ -> GERPpp
    -e '/^##INFO=<ID=GERP\+\+_NR([,>])/{s/ID=GERP\+\+_NR/ID=GERPpp_NR/;}' \
    -e 's/;GERP++_NR=/;GERPpp_NR=/g' \
    -e '/^##INFO=<ID=GERP\+\+_RS([,>])/{s/ID=GERP\+\+_RS/ID=GERPpp_RS/;}' \
    -e 's/;GERP++_RS=/;GERPpp_RS=/g' \
    "$vcf"
}

# Optional: add missing INFO headers for tags present in body (prevents "not defined in header" warnings)
add_missing_info_headers_if_needed() {
  local vcf="$1"
  local tmp="${vcf}.tmp"
  local hdr
  hdr="$(mktemp)"

  # Add only minimal, safe definitions.
  cat > "$hdr" <<'EOF'
##INFO=<ID=Alt,Number=1,Type=String,Description="Alt field from upstream annotation">
##INFO=<ID=Start,Number=1,Type=Integer,Description="Start field from upstream annotation">
EOF

  bcftools annotate -h "$hdr" -Ov -o "$tmp" "$vcf"
  mv "$tmp" "$vcf"
  rm -f "$hdr"
}

main() {
  local outdir
  outdir="$(dirname "${outputPrefix}")"
  assert_writable_dir "$outdir"

  # Output paths
  norm0Vcf="${outputPrefix}.biallelic.leftnorm0.vcf.gz"
  filteredVcf="${outputPrefix}.biallelic.leftnorm.vcf.gz"
  sitesVcf="${outputPrefix}.sites.vcf"

  annovarOutPrefix="${outputPrefix}.annovar"
  annovarMultiVcf="${annovarOutPrefix}.${refBuild}_multianno.vcf"
  annovarMultiVcfgz="${annovarMultiVcf}.gz"

  vcfAnnotated="${outputPrefix}.annovar.GT.vcf"
  finalVcfgz="${vcfAnnotated}.gz"

  # Clean stale outputs that can cause confusing permission/old-file reuse issues
  rm -f "${annovarOutPrefix}."* 2>/dev/null || true
  rm -f "${vcfAnnotated}" "${finalVcfgz}" "${finalVcfgz}.tbi" 2>/dev/null || true

  # 1) Normalize: split multiallelic -> biallelic, left-normalize
  bcftools norm -m-both -f "${refFASTA}" -Oz -o "${norm0Vcf}" "${vcfFile}"
  tabix -p vcf -f "${norm0Vcf}"

  # 2) Filter: KEEP only variants within BED regions (targets)
  bcftools filter -R "${bed}" -Oz -o "${filteredVcf}" "${norm0Vcf}"
  tabix -p vcf -f "${filteredVcf}"

  # 3) Sites-only VCF for ANNOVAR (keeps INFO, drops samples/FORMAT)
  bcftools view -G -Ov "${filteredVcf}" > "${sitesVcf}"

  # 4) ANNOVAR
  "${annovarFolder}/table_annovar.pl" "${sitesVcf}" "${annovardbFolder}" \
    -thread "${thread}" -buildver "${refBuild}" -out "${annovarOutPrefix}" \
    -remove -protocol "${protocol}" -operation "${operation}" \
    -nastring . -vcfinput > "${annovarOutPrefix}.log"

  # 5) Compress + index ANNOVAR multianno VCF
  bgzip -f "${annovarMultiVcf}"
  tabix -p vcf -f "${annovarMultiVcfgz}"

  # 6) Merge ANNOVAR INFO back into original multi-sample VCF (preserve GT/samples)
  bcftools annotate \
    -a "${annovarMultiVcfgz}" \
    -c INFO \
    --force \
    -Ov -o "${vcfAnnotated}" \
    "${filteredVcf}"

  # 7) Fix illegal INFO tag names (header + body keys)
  fix_vcf_info_tags_inplace "${vcfAnnotated}"

  # 8) Make INFO Types consistent (header-only edits)
  sed_inplace \
    -e '/^##INFO=<ID=gnomAD/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=ExAC_nontcga/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=CADD/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=REVEL/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=DANN_score/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=GME_AF/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=abraom_freq/s/Type=String/Type=Float/' \
    "${vcfAnnotated}"

  # 9) (Optional) silence "Alt/Start not defined in header" by adding header definitions
  add_missing_info_headers_if_needed "${vcfAnnotated}"

  # 10) Compress + index final output
  bgzip -f "${vcfAnnotated}"
  tabix -p vcf -f "${finalVcfgz}"

  # 11) Cleanup intermediates
  rm -f \
    "${sitesVcf}" \
    "${norm0Vcf}" \
    "${norm0Vcf}.tbi" \
    "${filteredVcf}.tbi" \
    "${annovarMultiVcfgz}" \
    "${annovarMultiVcfgz}.tbi" \
    "${annovarOutPrefix}.${refBuild}_multianno.txt" \
    "${annovarOutPrefix}.avinput" \
    "${annovarOutPrefix}.log"

  echo "Done."
  echo "Final: ${finalVcfgz}"
}

main
