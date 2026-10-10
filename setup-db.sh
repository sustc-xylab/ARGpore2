#!/usr/bin/env bash
# 只准备数据库：从原出处下载钉死的快照，不解压后删除源包。
# 不替换为 PLSDB 2025 等后续版本。403/404/文件过小视为错误，并说明手动放置路径。
# 依赖 setup-tools.sh 已经把 seqkit/blast/perl/lastdb 准备好。
set -euo pipefail

SCRIPT=$(realpath "$0")
DIR=$(dirname "$SCRIPT")
DB="$DIR/database"
NPROC="${NPROC:-$(nproc 2>/dev/null || echo 8)}"
CONDA_ENV="${ARGPORE_CONDA_ENV:-argpore2}"
FAILED=0

# 版本锁定（与 argpore.sh 及 2022 年文献所用参考数据一致，不是各数据库的最新发布）
PLSDB_URL="https://ndownloader.figshare.com/files/21961095"
PLSDB_PAGE="https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download"
CEN_URL="https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz"

log() { echo "[$(date +'%F %T')] $*"; }
warn() { echo "[WARNING] $*" >&2; }
err() { echo "[ERROR] $*" >&2; FAILED=1; }

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

if CONDA=$(find_conda); then
  if "$CONDA" env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
    CONDA_PREFIX="$("$CONDA" info --base)/envs/${CONDA_ENV}"
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
    unset PERL5LIB
  fi
fi

LASTDB="$DIR/bin/last-983/src/lastdb"
PERL="${CONDA_PREFIX:-}/bin/perl"
if [ ! -x "${PERL:-}" ]; then
  PERL=$(command -v perl || true)
fi
[ -x "$LASTDB" ] || { echo "[ERROR] 找不到 lastdb，请先: bash \"$DIR/setup-tools.sh\"" >&2; exit 1; }
[ -n "$PERL" ] && [ -x "$PERL" ] || { echo "[ERROR] 找不到 perl，请先: bash \"$DIR/setup-tools.sh\"" >&2; exit 1; }
if ! "$PERL" -MBio::SeqIO -e 1; then
  echo "[ERROR] perl 无 Bio::SeqIO，请先: bash \"$DIR/setup-tools.sh\"" >&2
  exit 1
fi
command -v wget >/dev/null 2>&1 || { echo "[ERROR] 需要 wget" >&2; exit 1; }
command -v curl >/dev/null 2>&1 || { echo "[ERROR] 需要 curl（用于探测 HTTP 状态）" >&2; exit 1; }

file_ok() {
  local f="$1" min="$2"
  [ -s "$f" ] && [ "$(stat -c%s "$f" 2>/dev/null || echo 0)" -ge "$min" ]
}

http_code() {
  curl -sI -L --max-time 30 -A "Mozilla/5.0" -o /dev/null -w "%{http_code}" "$1" || echo "000"
}

# 可续传；失败返回 1。不把 403 当成成功。
download() {
  local url="$1" out="$2"
  log "下载 $url -> $out"
  if wget -c --tries=20 --retry-connrefused --timeout=60 --waitretry=10 \
      --progress=bar:force:noscroll \
      -O "$out" "$url"; then
    return 0
  fi
  warn "wget 失败: $url"
  return 1
}

need_index() {
  local prefix="$1"
  [ ! -e "${prefix}.prj" ] && [ ! -e "${prefix}.bck" ]
}

mkdir -p "$DB"
cd "$DB"

# --- markers（git LFS，仓内源包） ---
if command -v git-lfs >/dev/null 2>&1; then
  log "git lfs pull（markers 等仓内大文件）"
  git -C "$DIR" lfs install
  git -C "$DIR" lfs pull || err "git lfs pull 失败。可稍后重试；markers 源包: $DB/markers.fasta.tar.xz"
else
  err "缺少 git-lfs。conda install -n $CONDA_ENV -c conda-forge git-lfs 后重跑 setup-db.sh"
fi

if [ ! -s markers.fasta ]; then
  if [ -s markers.fasta.tar.xz ]; then
    log "解压 markers.fasta.tar.xz（保留源包）"
    tar -jxvf markers.fasta.tar.xz
  else
    err "缺少 markers.fasta.tar.xz。请确认 clone 后执行了 git lfs pull"
  fi
fi
if [ -s markers.fasta ] && need_index markers.lastindex; then
  log "lastdb markers (-P $NPROC)"
  "$LASTDB" -Q 0 markers.lastindex markers.fasta -P "$NPROC"
fi
if [ -s markers.fasta ] && [ ! -s markers.fasta.length ]; then
  "$PERL" "$DIR/bin/fastaNameLengh.pl" markers.fasta > markers.fasta.length
fi

# --- SARG（仓内源包，版本 SARG_20211207） ---
if [ ! -s SARG_20211207_14210_filter.ffn ]; then
  if [ -s SARG_20211207_14210_filter.ffn.tar.gz ]; then
    log "解压 SARG（保留源包）"
    tar -zxvf SARG_20211207_14210_filter.ffn.tar.gz
  else
    err "缺少 SARG_20211207_14210_filter.ffn.tar.gz（应在 git clone 的 database/ 中）"
  fi
fi
if [ -s SARG_20211207_14210_filter.ffn ] && need_index SARG_20211207_14210_filter.ffn; then
  log "lastdb SARG (-P $NPROC)"
  "$LASTDB" -Q 0 SARG_20211207_14210_filter.ffn SARG_20211207_14210_filter.ffn -P "$NPROC"
fi
if [ -s SARG_20211207_14210_filter.ffn ] && [ ! -s SARG_20211207_14210_filter.ffn.length ]; then
  "$PERL" "$DIR/bin/fastaNameLengh.pl" SARG_20211207_14210_filter.ffn > SARG_20211207_14210_filter.ffn.length
fi

# --- ESCG ---
if [ ! -s ESCG.fna ]; then
  if [ -s ESCG.fna.tar.xz ]; then
    log "解压 ESCG（保留源包）"
    tar -jxvf ESCG.fna.tar.xz
  else
    err "缺少 ESCG.fna.tar.xz"
  fi
fi
if [ -s ESCG.fna ] && need_index ESCG.fna; then
  log "lastdb ESCG (-P $NPROC)"
  "$LASTDB" -Q 0 ESCG.fna ESCG.fna -P "$NPROC"
fi
if [ -s ESCG.fna ] && [ ! -s ESCG.fna.length ]; then
  "$PERL" "$DIR/bin/fastaNameLengh.pl" ESCG.fna > ESCG.fna.length
fi

# --- lineage 2020-06-16 ---
if [ ! -s 2020-06-16_lineage.tab ]; then
  if [ -s 2020-06-16_lineage.tab.tar.xz ]; then
    log "解压 lineage（保留源包）"
    tar -jxvf 2020-06-16_lineage.tab.tar.xz
    if [ -f database/2020-06-16_lineage.tab ]; then
      mv database/2020-06-16_lineage.tab .
      rm -rf database
    fi
  else
    err "缺少 2020-06-16_lineage.tab.tar.xz"
  fi
fi

# --- PLSDB 2020_03_04（Figshare 固定文件，禁止改下 2025） ---
if file_ok PLSDB_2020_03_04.fna 500000000; then
  log "已有 PLSDB_2020_03_04.fna，跳过下载"
elif [ -d PLSDB_2020_03_04 ]; then
  log "从已有 BLAST 库导出 FASTA"
  "$DIR/bin/ncbi-blast+/bin/blastdbcmd" \
    -db PLSDB_2020_03_04/plsdb.fna -dbtype nucl -entry all -outfmt "%f" \
    -out PLSDB_2020_03_04.fna || err "blastdbcmd 导出 PLSDB FASTA 失败"
else
  code=$(http_code "$PLSDB_URL")
  log "探测 Figshare PLSDB HTTP $code"
  if [ "$code" = "403" ] || [ "$code" = "404" ] || [ "$code" = "000" ]; then
    err "无法从 Figshare 下载 PLSDB 2020_03_04（HTTP $code）。
  所需文件为该日期的发布版本，不会替换为 PLSDB 2025。
  官方历史版本页: $PLSDB_PAGE （列表中的 version 2020_03_04）
  Figshare 文件: $PLSDB_URL
  请在可访问 Figshare 的网络下载后，把文件放到:
    $DB/PLSDB.zip
  或已解压的 BLAST 库目录:
    $DB/PLSDB_2020_03_04/
  或现成 FASTA（仍用此文件名，argpore.sh 写死了）:
    $DB/PLSDB_2020_03_04.fna
  然后重新运行: bash \"$DIR/setup-db.sh\""
  else
    if [ -s PLSDB.zip ] && file_ok PLSDB.zip 100000000; then
      log "使用已有 PLSDB.zip"
    elif download "$PLSDB_URL" PLSDB.zip && file_ok PLSDB.zip 100000000; then
      log "已下载 PLSDB.zip（保留，不删除）"
    else
      err "PLSDB.zip 下载不完整。请手动放置 $DB/PLSDB.zip 或 $DB/PLSDB_2020_03_04.fna 后重跑。
  链接: $PLSDB_URL
  版本页: $PLSDB_PAGE"
    fi
    if file_ok PLSDB.zip 100000000 && [ ! -d PLSDB_2020_03_04 ]; then
      log "解压 PLSDB.zip"
      unzip -o PLSDB.zip
      if [ -d data/pls/PLSDB_2020_03_04 ]; then
        mv data/pls/PLSDB_2020_03_04 .
      fi
      rm -rf data
    fi
    if [ -d PLSDB_2020_03_04 ] && [ ! -s PLSDB_2020_03_04.fna ]; then
      "$DIR/bin/ncbi-blast+/bin/blastdbcmd" \
        -db PLSDB_2020_03_04/plsdb.fna -dbtype nucl -entry all -outfmt "%f" \
        -out PLSDB_2020_03_04.fna || err "blastdbcmd 导出 PLSDB FASTA 失败"
    fi
  fi
fi

if file_ok PLSDB_2020_03_04.fna 500000000; then
  if [ ! -s PLSDB_2020_03_04.fna.length ]; then
    grep ">" PLSDB_2020_03_04.fna | sed 's/>//' | sed 's/ /\t/' > PLSDB_2020_03_04.fna.name
    "$PERL" "$DIR/bin/fastaNameLengh.pl" PLSDB_2020_03_04.fna > PLSDB_2020_03_04.fna.length
  fi
  if need_index PLSDB_2020_03_04.fna.lastindex; then
    log "lastdb PLSDB（较慢，-P $NPROC）"
    "$LASTDB" -Q 0 PLSDB_2020_03_04.fna.lastindex PLSDB_2020_03_04.fna -P "$NPROC"
  else
    log "已有 PLSDB last 索引，跳过"
  fi
else
  if [ "$FAILED" -eq 0 ]; then
    err "缺少有效的 PLSDB_2020_03_04.fna（需约 2.1G）。请按上面的路径手动放置后重跑 setup-db.sh"
  fi
fi

# --- Centrifuge p_compressed+h+v（约 2016 官方索引，S3） ---
if file_ok p_compressed+h+v.1.cf 1000000000; then
  log "已有 centrifuge 索引，跳过下载"
else
  code=$(http_code "$CEN_URL")
  log "探测 Centrifuge 索引 HTTP $code"
  if [ "$code" = "403" ] || [ "$code" = "404" ] || [ "$code" = "000" ]; then
    err "无法下载 Centrifuge 索引（HTTP $code）。
  链接: $CEN_URL
  请将解压后的三个文件放到:
    $DB/p_compressed+h+v.1.cf
    $DB/p_compressed+h+v.2.cf
    $DB/p_compressed+h+v.3.cf
  然后重新运行 setup-db.sh。不完整的 p+b+v.tar.gz 请删掉后重下。"
  else
    if [ -s p+b+v.tar.gz ] && ! file_ok p+b+v.tar.gz 1000000000; then
      warn "已有 p+b+v.tar.gz 但体积过小，尝试续传"
    fi
    if download "$CEN_URL" p+b+v.tar.gz && file_ok p+b+v.tar.gz 1000000000; then
      log "解压 centrifuge 索引（保留 tar.gz）"
      tar -zxvf p+b+v.tar.gz
    else
      err "Centrifuge 索引下载不完整（常见于中断）。请 wget -c 续传或手动放置三个 .cf 后重跑。
  $CEN_URL -> $DB/p+b+v.tar.gz"
    fi
  fi
fi

if ! file_ok p_compressed+h+v.1.cf 1000000000; then
  if [ "$FAILED" -eq 0 ]; then
    err "缺少有效的 p_compressed+h+v.1.cf"
  fi
fi

echo
if [ "$FAILED" -ne 0 ]; then
  echo "[ERROR] 数据库未准备完整。已完成的步骤会保留，修好文件后可重跑 setup-db.sh。" >&2
  echo "检查命令: bash \"$DIR/check-setup.sh\"" >&2
  exit 1
fi

log "数据库准备完成（版本锁定为 2020–2021 年所用参考数据，不是各数据库的最新发布）"
bash "$DIR/check-setup.sh"
