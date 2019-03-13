#!/bin/sh

set -e

test $# -ge 1

cd "`dirname "$0"`"

name="$1"

make "$name"
./"$@"

export GCOV_PREFIX="`pwd`"
export GCOV_PREFIX_STRIP="`echo "$GCOV_PREFIX" | tr -Cd / | wc -c`"

mkdir -p "$name-data"
mv -f "$name.gc"* "$name-data"
lcov -t "$name" -c -d "$name-data" -o "$name.info"
rm -fr "$name-coverage"
genhtml -o "$name-coverage" "$name.info"
