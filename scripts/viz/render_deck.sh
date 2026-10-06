#!/bin/bash
# Render a .pptx to per-slide JPEGs for visual QA.
#
# LibreOffice is NOT in the system path on Woody and cannot be installed with
# dnf (no sudo). The official 26.x RPMs need glibc 2.34 and this host has 2.28,
# so 7.6.7.2 is used -- the last line that still builds against the RHEL 8
# baseline. Installed without root by extracting the RPMs with rpm2cpio into
# ~/software/libreoffice/root76.
#
# To reinstall from scratch:
#   V=7.6.7.2; D=~/software/libreoffice; mkdir -p $D/src76 && cd $D/src76
#   curl -sLO https://downloadarchive.documentfoundation.org/libreoffice/old/$V/rpm/x86_64/LibreOffice_${V}_Linux_x86-64_rpm.tar.gz
#   tar xzf LibreOffice_${V}_Linux_x86-64_rpm.tar.gz
#   mkdir -p $D/root76 && cd $D/root76
#   for r in $D/src76/*/RPMS/*.rpm; do rpm2cpio "$r" | cpio -idm --quiet; done
#
#   render_deck.sh deck.pptx [outdir]
set -euo pipefail
SOFF=${SOFF:-$HOME/software/libreoffice/root76/opt/libreoffice7.6/program/soffice}
[[ -x "$SOFF" ]] || { echo "FATAL: no soffice at $SOFF -- see the header" >&2; exit 1; }
DECK=$(readlink -f "$1"); OUT=$(readlink -f "${2:-$(dirname "$DECK")}")
mkdir -p "$OUT"; cd "$OUT"
rm -f slide-*.jpg
# -env:UserInstallation keeps it out of a real profile and lets it run headless
"$SOFF" --headless --norestore \
  -env:UserInstallation=file://"$HOME"/.config/libreoffice-headless \
  --convert-to pdf "$DECK" --outdir "$OUT" >/dev/null 2>&1
pdftoppm -jpeg -r "${DPI:-100}" "$OUT/$(basename "${DECK%.pptx}").pdf" slide
ls -1 "$OUT"/slide-*.jpg
