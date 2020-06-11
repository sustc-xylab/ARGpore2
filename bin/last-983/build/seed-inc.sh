#! /bin/sh

# This generates source code from subset seeds.

cat <<EOF
const struct {
  const char *name;
  const char *text;
} subsetSeeds[] = {
EOF
for i in "$@"
do
    basename $i .seed | sed 's/.*/{"&", "\\/'
    grep -v '^# ' $i | awk NF | sed 's/$/\\n\\/'
    echo '"},'
    echo
done
echo "};"
