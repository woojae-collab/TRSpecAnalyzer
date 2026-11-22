#!/bin/bash
# TRSpecAnalyzer, Python version - installation on macOS.
# Double-click this file (or, in Terminal:  bash install_mac.command).
# It makes the folder .venv inside this folder and downloads numpy,
# scipy, matplotlib and PySide6 into it. Nothing outside this folder is
# changed.

cd "$(dirname "$0")" || exit 1

finish() {          # $1 = exit code
    echo
    read -n 1 -s -r -p "Press any key to finish."
    echo
    exit "$1"
}

fail() {
    echo
    echo "============================================================"
    echo " Installation did not finish. Please copy the text above"
    echo " and send it along."
    echo "============================================================"
    finish 1
}

usable() {          # is $1 a Python 3.10 or later, 64-bit?
    [ -x "$1" ] || return 1
    "$1" -c "import sys; sys.exit(0 if sys.version_info[:2] >= (3, 10) and sys.maxsize > 2**32 else 1)" >/dev/null 2>&1
}

echo "============================================================"
echo " TRSpecAnalyzer, Python version - installation"
echo "============================================================"
echo

# Where a Python can be: the installer of python.org, Homebrew, the PATH.
# (/usr/bin/python3 is left out: without Apple's developer tools it only
# opens a dialog that offers to install them, and with them it is too
# old.)
PY=""
for v in 3.13 3.12 3.14 3.11 3.10; do
    c="/Library/Frameworks/Python.framework/Versions/$v/bin/python3"
    if usable "$c"; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
    for c in /opt/homebrew/bin/python3 /usr/local/bin/python3; do
        if usable "$c"; then PY="$c"; break; fi
    done
fi
if [ -z "$PY" ]; then
    c="$(command -v python3 2>/dev/null)"
    if [ -n "$c" ] && [ "$c" != "/usr/bin/python3" ] && usable "$c"; then PY="$c"; fi
fi

if [ -z "$PY" ]; then
    echo "A Python that can be used (3.10 or later) was not found."
    echo
    echo "1. Download Python 3.13 for macOS from"
    echo "   https://www.python.org/downloads/  and install it."
    echo "2. Open this file again."
    finish 1
fi

echo "Python found:"
"$PY" --version
echo "  ($PY)"

if [ ! -x ".venv/bin/python" ]; then
    echo
    echo "Creating a private environment in the folder .venv ..."
    "$PY" -m venv .venv || fail
fi

echo
echo "Installing numpy, scipy, matplotlib and PySide6."
echo "This needs an internet connection and takes a few minutes ..."
echo
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt || fail
.venv/bin/python -c "import numpy, scipy, matplotlib, PySide6; print('Packages are in place.')" || fail

echo
echo "============================================================"
echo " Installation finished."
echo " Start the program with run_mac.command"
echo "============================================================"
finish 0
