#!/bin/sh

cd "$( dirname "${BASH_SOURCE[0]}")"

gzip -c S288cvsRM11_1a.maf > S288cvsRM11_1a.maf.gz

../multiz M=10 S288cvsRM11_1a.maf \
    S288cvsSpar.maf \
    1 out1 out2 \
    | grep -v "#" > maf.maf

../multiz M=10 S288cvsRM11_1a.maf.gz \
    S288cvsSpar.maf \
    1 out1 out2 \
    | grep -v "#" > maf.gz.maf

cmp --silent maf.maf maf.gz.maf || echo "files are different"

rm out1 out2
