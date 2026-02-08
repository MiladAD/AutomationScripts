#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: pipeline.sh <vcfFile> <outputPrefix> <refFASTA> <bed> <annovarFolder> <annovardbFolder> <refBuild> <thread> [protocol] [operation]

  vcfFile:         Input VCF file (optionally bgzipped)
  outputPrefix:    Prefix for output files (can include output dir)
  refFASTA:        Reference FASTA file (indexed for bcftools norm)
  bed:             BED file (regions to KEEP with bcftools -R)
  annovarFolder:   Path to ANNOVAR scripts (contains table_annovar.pl)
  annovardbFolder: Path to ANNOVAR databases
  refBuild:        Reference build (e.g., hg19 or hg38)
  thread:          ANNOVAR threads
  protocol:        (Optional) ANNOVAR protocol string
  operation:       (Optional) ANNOVAR operation string

Notes:
  - The BED is used as "target regions to keep" (bcftools filter -R).
  - ANNOVAR runs on a sites-only VCF, then its INFO annotations are merged back
    into the original multi-sample VCF to preserve GT/sample columns.
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

# Default protocol/operation (hg19 set by default; you can override via args 9/10)
protocol="refGeneWithVer,knownGene,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad211_genome,gnomad211_exome,exac03,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
operation="g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"

# If you want hg38 defaults, pass custom protocol/operation OR uncomment these and comment hg19 above:
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

# Handle sed -i portability (GNU vs BSD/macOS)
sed_inplace() {
  if sed --version >/dev/null 2>&1; then
    # GNU sed
    sed -i "$@"
  else
    # BSD sed (macOS)
    sed -i '' "$@"
  fi
}

main() {
  # Outputs
  norm0Vcf="${outputPrefix}.biallelic.leftnorm0.vcf.gz"
  filteredVcf="${outputPrefix}.biallelic.leftnorm.vcf.gz"
  sitesVcf="${outputPrefix}.sites.vcf"
  annovarOutPrefix="${outputPrefix}.annovar"
  annovarMultiVcf="${annovarOutPrefix}.${refBuild}_multianno.vcf"
  annovarMultiVcfgz="${annovarMultiVcf}.gz"
  vcfAnnotated="${outputPrefix}.annovar.GT.vcf"
  finalVcfgz="${vcfAnnotated}.gz"

  # 1) Split multiallelic -> biallelic, left-normalize
  bcftools norm -m-both -f "${refFASTA}" -Oz -o "${norm0Vcf}" "${vcfFile}"
  tabix -p vcf -f "${norm0Vcf}"

  # 2) Keep only variants within BED regions (targets)
  bcftools filter -R "${bed}" -Oz -o "${filteredVcf}" "${norm0Vcf}"
  tabix -p vcf -f "${filteredVcf}"

  # 3) Create a sites-only VCF for ANNOVAR (removes samples/FORMAT to speed up)
  bcftools view -G -Ov "${filteredVcf}" > "${sitesVcf}"

  # 4) Run ANNOVAR
  "${annovarFolder}/table_annovar.pl" "${sitesVcf}" "${annovardbFolder}" \
    -thread "${thread}" -buildver "${refBuild}" -out "${annovarOutPrefix}" \
    -remove -protocol "${protocol}" -operation "${operation}" \
    -nastring . -vcfinput > "${annovarOutPrefix}.log"

  # 5) Compress + index ANNOVAR multianno VCF
  bgzip -f "${annovarMultiVcf}"
  tabix -p vcf -f "${annovarMultiVcfgz}"

  # 6) Merge ANNOVAR INFO back into the ORIGINAL multi-sample VCF (preserve GT/samples)
  bcftools annotate \
    -a "${annovarMultiVcfgz}" \
    -c INFO \
    --force \
    -Ov -o "${vcfAnnotated}" \
    "${filteredVcf}"

  # 7) Header/body string normalizations (keep Type edits limited to header lines)
  #    WARNING: The global substitutions below affect body too; kept because you had them,
  #    but they are intentionally targeted at known problematic tokens.
  sed_inplace \
    -e '/^##INFO=<ID=gnomAD/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=ExAC_nontcga/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=CADD/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=REVEL/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=DANN_score/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=GME_AF/s/Type=String/Type=Float/' \
    -e '/^##INFO=<ID=abraom_freq/s/Type=String/Type=Float/' \
    -e 's/Eigen-/Eigen_/g' \
    -e 's/GERP++/GERPpp/g' \
    -e 's/PC-/PC_/g' \
    -e 's/M-CAP/M_CAP/g' \
    -e 's/fathmm-/fathmm_/g' \
    "${vcfAnnotated}"

  # 8) Compress + index final VCF
  bgzip -f "${vcfAnnotated}"
  tabix -p vcf -f "${finalVcfgz}"

  # 9) Cleanup intermediates
  rm -f \
    "${sitesVcf}" \
    "${norm0Vcf}" \
    "${norm0Vcf}.tbi" \
    "${filteredVcf}.tbi" \
    "${annovarMultiVcfgz}" \
    "${annovarMultiVcfgz}.tbi" \
    "${annovarOutPrefix}.${refBuild}_multianno.txt" \
    "${annovarOutPrefix}.avinput"

  echo "Done."
  echo "Final: ${finalVcfgz}"
}

main
