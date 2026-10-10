#!/usr/bin/env bash
# 只装/链接运行时工具，不下数据库。数据库见 setup-db.sh。
# argpore.sh 仍读取 $DIR/bin/seqkit、bin/ncbi-blast+、bin/centrifuge。
set -euo pipefail

SCRIPT=$(realpath "$0")
DIR=$(dirname "$SCRIPT")
CONDA_ENV="${ARGPORE_CONDA_ENV:-argpore2}"
R_ENV="${ARGPORE_R_ENV:-argpore2-r}"

log() { echo "[$(date +'%F %T')] $*"; }
warn() { echo "[WARNING] $*" >&2; }
die() { echo "[ERROR] $*" >&2; exit 1; }

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
    if [ -x "$c" ]; then
      echo "$c"
      return 0
    fi
  done
  return 1
}

CONDA=$(find_conda) || die "找不到 conda。请先安装 miniforge/miniconda，并把 conda 加入 PATH。"

if ! "$CONDA" env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
  die "conda 环境 '$CONDA_ENV' 不存在。请先创建：
  conda env create -f \"$DIR/environment-tools.yml\"
  # 或: conda create -n $CONDA_ENV -y -c conda-forge -c bioconda seqkit blast centrifuge parallel perl-bioperl git-lfs wget"
fi

CONDA_PREFIX="$("$CONDA" info --base)/envs/${CONDA_ENV}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"
unset PERL5LIB

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "缺少命令: $1（应在 conda env $CONDA_ENV 中）"
}

need_cmd seqkit
need_cmd blastdbcmd
need_cmd blastn
need_cmd centrifuge
need_cmd parallel
need_cmd wget

PERL="${CONDA_PREFIX}/bin/perl"
[ -x "$PERL" ] || die "找不到 $PERL"
if ! "$PERL" -MBio::SeqIO -e 1; then
  die "Bio::SeqIO 不可用。请: conda install -n $CONDA_ENV -c bioconda perl-bioperl"
fi

LASTDB="$DIR/bin/last-983/src/lastdb"
[ -x "$LASTDB" ] || die "找不到 lastdb: $LASTDB（仓库自带 last-983，不应缺失）"

PREFIX="$(dirname "$(command -v seqkit)")"
log "把 conda 工具软链到 $DIR/bin（argpore.sh 的硬编码路径）"
ln -sfn "$PREFIX/seqkit" "$DIR/bin/seqkit"
mkdir -p "$DIR/bin/ncbi-blast+/bin"
ln -sfn "$PREFIX/blastdbcmd" "$DIR/bin/ncbi-blast+/bin/blastdbcmd"
ln -sfn "$PREFIX/blastn" "$DIR/bin/ncbi-blast+/bin/blastn"
mkdir -p "$DIR/bin/centrifuge"
ln -sfn "$PREFIX/centrifuge" "$DIR/bin/centrifuge/centrifuge"
if command -v centrifuge-build >/dev/null 2>&1; then
  ln -sfn "$PREFIX/centrifuge-build" "$DIR/bin/centrifuge/centrifuge-build"
fi
if command -v centrifuge-inspect >/dev/null 2>&1; then
  ln -sfn "$PREFIX/centrifuge-inspect" "$DIR/bin/centrifuge/centrifuge-inspect"
fi

"$DIR/bin/seqkit" version
"$DIR/bin/ncbi-blast+/bin/blastdbcmd" -version | head -n 1 || true
"$DIR/bin/centrifuge/centrifuge" --version | head -n 1 || true
log "lastdb: $LASTDB"
log "perl Bio::SeqIO: OK ($PERL)"

if command -v git-lfs >/dev/null 2>&1; then
  log "git-lfs: $(command -v git-lfs)"
else
  warn "未找到 git-lfs。setup-db.sh 拉取 markers 时需要它: conda install -n $CONDA_ENV -c conda-forge git-lfs"
fi

if "$CONDA" env list | awk '{print $1}' | grep -qx "$R_ENV"; then
  if "$CONDA" run -n "$R_ENV" Rscript -e 'for (p in c("plyr","data.table","dplyr","foreach","doParallel")) library(p, character.only=TRUE)' >/dev/null 2>&1; then
    log "R env $R_ENV: packages OK"
  else
    warn "conda 环境 $R_ENV 存在，但缺少 plyr/data.table/dplyr/foreach/doParallel"
  fi
else
  warn "未找到 R 环境 '$R_ENV'。运行 argpore.sh 需要 R。请: conda env create -f \"$DIR/environment-r.yml\""
  if command -v Rscript >/dev/null 2>&1; then
    log "当前 PATH 中有 Rscript: $(command -v Rscript)"
  fi
fi

log "工具安装完成。下一步: bash \"$DIR/setup-db.sh\""
