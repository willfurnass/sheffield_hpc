#!/bin/bash

SRC_FILE=newfloat.dtx
OUTPUT_FILE=newfloat.sty
SRC_FILE_URL=http://mirrors.ctan.org/macros/latex/contrib/caption/$SRC_FILE

function usage {
    echo 1>&2 "Download $OUTPUT_FILE and save to a local directory"
    echo 1>&2 "usage: $0 output_dirname"
    exit 1
}

[[ $# -eq 1 ]] || usage
[[ -d "$1" ]] || usage
OUTPUT_DIR="$1"

mkdir .latexbuild
pushd .latexbuild
wget -N "$SRC_FILE_URL"

cat << EOF > buildnewfloat.ins
\def\batchfile{buildnewfloat.ins}
\input docstrip
\generate{\file{${OUTPUT_FILE}}{\from{${SRC_FILE}}{package}}}
EOF

latex buildnewfloat.ins
popd
cp ".latexbuild/$OUTPUT_FILE" "$OUTPUT_DIR"
echo "$OUTPUT_FILE saved to $OUTPUT_DIR"
rm -r .latexbuild
