#!/bin/bash
# TRSpecAnalyzer, Python version - start on macOS.
# Double-click this file (or, in Terminal:  bash run_mac.command).

cd "$(dirname "$0")" || exit 1

if [ ! -x ".venv/bin/python" ]; then
    echo "The program is not installed yet. Open install_mac.command first."
    echo
    read -n 1 -s -r -p "Press any key to finish."
    echo
    exit 1
fi

echo "TRSpecAnalyzer, Python version. This window shows messages of the"
echo "program. It can be closed after the program is closed."
.venv/bin/python -m trspecgui "$@"
rc=$?
if [ "$rc" -ne 0 ]; then
    echo
    echo "The program ended with an error. Please copy the text above and"
    echo "send it along."
    echo
    read -n 1 -s -r -p "Press any key to finish."
    echo
    exit 1
fi
exit 0
