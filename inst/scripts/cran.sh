#!/usr/bin/env bash
#SBATCH --job-name=_cran
#SBATCH --account=PETERS-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --mem=28800
#SBATCH --time=12:00:00
#SBATCH --error=/home/jhz22/R/work/_cran_%A_%a.err
#SBATCH --output=/home/jhz22/R/work/_cran_%A_%a.out
#SBATCH --export=ALL

set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl ceuadmin/R mono/5.0.1.1 texlive

export TMPDIR=/rds/user/jhz22/hpc-work/work

src="$HOME/pQTLtools"
dst="$HOME/R/pQTLtools"
log="$HOME/work/pQTLtools_cran.log"

exec > >(tee -a "$log") 2>&1
echo "=== CRAN pipeline start $(date) on $(hostname) ==="

rm -rf "$dst"
mkdir -p "$dst"
rsync -a --delete \
  --exclude='docs/' --exclude='pkgdown/' \
  --exclude='README.Rmd' \
  --exclude='.*' "$src/" "$dst/"

cd "$HOME"

ver=$(awk '/^Version:/ {print $2}' "$src/DESCRIPTION")
pkg="pQTLtools_${ver}.tar.gz"
[ -n "${ver:-}" ] || { echo "Version not found"; exit 1; }

R CMD build --resave-data --compact-vignettes=both "$src"
R CMD INSTALL "$pkg"
R CMD check --as-cran --run-donttest "$pkg"

rm -rf "$HOME/R/pQTLtools.Rcheck" || true
mv "$pkg" "$HOME/R/"
mv pQTLtools.Rcheck "$HOME/R/"

echo "=== CRAN pipeline finished $(date) ==="
