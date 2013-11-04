#!/bin/bash

# note:
# - drivcav set: download 11M, uncompressed 35M
# - fidap set:   download  1M, uncompressed  5M

for URL in \
  ftp://math.nist.gov/pub/MatrixMarket2/SPARSKIT/drivcav/{index.html,e{05,40}r{0100,0500}{.html,{,_rhs1}.mtx.gz}} \
  ftp://math.nist.gov/pub/MatrixMarket2/SPARSKIT/fidap/{index.html,fidapm{03,13,33}{.html,_info.txt,{,_rhs1}.mtx.gz}} \
  ; do
  wget -nv -nH -x $URL
done
mv -f pub/MatrixMarket2/SPARSKIT/* ./ && rm -rf pub

gunzip -f */*.gz

