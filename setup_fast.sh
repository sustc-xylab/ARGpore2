#!/usr/bin/env bash
# 兼容入口：改为依次调用 setup-tools.sh 与 setup-db.sh。
set -euo pipefail
DIR=$(dirname "$(realpath "$0")")
echo "setup_fast.sh 已拆成 setup-tools.sh + setup-db.sh，正在转调..."
bash "$DIR/setup-tools.sh"
bash "$DIR/setup-db.sh"
