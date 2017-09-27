#!/bin/sh
root -l -b <<EOF
.L myHists.cxx+
.L e_AnalyzeTree.cxx
.x e_AnaTree.C("$1")
.q
EOF
