#!/usr/bin/env bash
# 检查工具软链与钉死数据库是否齐，缺则非 0 退出。
# 用法: bash check-setup.sh
#       bash argpore.sh --check
set -euo pipefail

SCRIPT=$(realpath "$0")
DIR=$(dirname "$SCRIPT")
DB="$DIR/database"
CONDA_ENV="${ARGPORE_CONDA_ENV:-argpore2}"
R_ENV="${ARGPORE_R_ENV:-argpore2-r}"
FAILED=0

ok() { echo "  OK    $*"; }
bad() { echo "  MISS  $*" >&2; FAILED=1; }
warn() { echo "  WARN  $*" >&2; }

find_conda() {
  if command -v conda >/dev/null 2>&1; then
    command -v conda
    return 0
  fi
  local c
  for c in \
    "${HOME}/miniforge3/condabin/conda" \
    "${HOME}/mambaforge/condabin/conda" \
    "${HOME}/miniconda3/condabin/conda" \
    "${HOME}/anaconda3/condabin/conda"
  do
    [ -x "$c" ] && echo "$c" && return 0
  done
  return 1
}

file_ok() {
  local f="$1" min="${2:-1}"
  [ -e "$f" ] && [ "$(stat -c%s "$f" 2>/dev/null || echo 0)" -ge "$min" ]
}

echo "ARGpore2 check  $DIR"
echo "--- tools ---"

if [ -x "$DIR/bin/seqkit" ] || [ -L "$DIR/bin/seqkit" ]; then
  ok "seqkit  $DIR/bin/seqkit"
else
  bad "seqkit  请运行 bash setup-tools.sh"
fi
if [ -x "$DIR/bin/ncbi-blast+/bin/blastdbcmd" ]; then
  ok "blastdbcmd"
else
  bad "blastdbcmd  $DIR/bin/ncbi-blast+/bin/blastdbcmd"
fi
if [ -x "$DIR/bin/centrifuge/centrifuge" ]; then
  ok "centrifuge"
else
  bad "centrifuge  $DIR/bin/centrifuge/centrifuge"
fi
if [ -x "$DIR/bin/last-983/src/lastdb" ]; then
  ok "lastdb"
else
  bad "lastdb  $DIR/bin/last-983/src/lastdb"
fi
if command -v parallel >/dev/null 2>&1; then
  ok "parallel  $(command -v parallel)"
else
  bad "GNU parallel 不在 PATH"
fi

if CONDA=$(find_conda 2>/dev/null); then
  CONDA_PREFIX="$("$CONDA" info --base)/envs/${CONDA_ENV}"
  if [ -x "${CONDA_PREFIX}/bin/perl" ] && "${CONDA_PREFIX}/bin/perl" -MBio::SeqIO -e 1 2>/dev/null; then
    ok "perl Bio::SeqIO  ($CONDA_ENV)"
  elif perl -MBio::SeqIO -e 1 2>/dev/null; then
    ok "perl Bio::SeqIO  (PATH)"
  else
    bad "perl Bio::SeqIO"
  fi
  if "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$R_ENV"; then
    if "$CONDA" run -n "$R_ENV" Rscript -e 'for (p in c("plyr","data.table","dplyr","foreach","doParallel")) library(p, character.only=TRUE)' >/dev/null 2>&1; then
      ok "R packages  ($R_ENV)"
    else
      bad "R 环境 $R_ENV 缺少 plyr/data.table/dplyr/foreach/doParallel"
    fi
  elif command -v Rscript >/dev/null 2>&1; then
    if Rscript -e 'for (p in c("plyr","data.table","dplyr","foreach","doParallel")) library(p, character.only=TRUE)' >/dev/null 2>&1; then
      ok "R packages  (PATH Rscript)"
    else
      warn "Rscript 在 PATH 中，但缺包；建议 conda env create -f environment-r.yml"
    fi
  else
    bad "Rscript  请 conda env create -f \"$DIR/environment-r.yml\""
  fi
else
  warn "找不到 conda，跳过 env 内 perl/R 检查"
  perl -MBio::SeqIO -e 1 2>/dev/null && ok "perl Bio::SeqIO  (PATH)" || bad "perl Bio::SeqIO"
fi

echo "--- databases (pinned snapshots, not latest) ---"
check_file() {
  local label="$1" path="$2" min="${3:-1}"
  if file_ok "$path" "$min"; then
    ok "$label"
  else
    bad "$label  $path"
  fi
}
check_index() {
  local label="$1" prefix="$2"
  if file_ok "${prefix}.prj" 1 || file_ok "${prefix}.bck" 1; then
    ok "$label"
  else
    bad "$label  ${prefix}.prj"
  fi
}
check_file "markers.fasta" "$DB/markers.fasta" 1000
check_index "markers lastindex" "$DB/markers.lastindex"
check_file "SARG_20211207" "$DB/SARG_20211207_14210_filter.ffn" 1000
check_index "SARG lastindex" "$DB/SARG_20211207_14210_filter.ffn"
check_file "ESCG.fna" "$DB/ESCG.fna" 1000
check_index "ESCG lastindex" "$DB/ESCG.fna"
check_file "lineage 2020-06-16" "$DB/2020-06-16_lineage.tab" 1000
check_file "structure.RData" "$DB/structure.RData" 1
check_file "taxa.info.RData" "$DB/taxa.info.RData" 1
check_file "PLSDB_2020_03_04.fna (not 2025)" "$DB/PLSDB_2020_03_04.fna" 500000000
check_index "PLSDB lastindex" "$DB/PLSDB_2020_03_04.fna.lastindex"
check_file "centrifuge p_compressed+h+v.1.cf" "$DB/p_compressed+h+v.1.cf" 1000000000
check_file "centrifuge .2.cf" "$DB/p_compressed+h+v.2.cf" 1000
check_file "centrifuge .3.cf" "$DB/p_compressed+h+v.3.cf" 1

echo
if [ "$FAILED" -ne 0 ]; then
  echo "check FAILED。工具: bash setup-tools.sh ；数据库: bash setup-db.sh" >&2
  exit 1
fi
echo "check OK"
exit 0
