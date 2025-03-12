#!/usr/bin/env bash

set -e
src="$HOME/pQTLtools"
dst="$HOME/R/pQTLtools"
log_file="$HOME/work/pQTLtools_copy.log"

remove_destination() {
  if [ -d "$dst" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Removing existing destination directory: $dst" | tee -a "$log_file"
    rm -rf "$dst"
  elif [ -e "$dst" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Removing existing file at destination: $dst" | tee -a "$log_file"
    rm -f "$dst"
  fi
}
{
  echo "----------------------------------------"
  echo "Copy Operation Started: $(date '+%Y-%m-%d %H:%M:%S')"
  remove_destination
  echo "Creating destination directory: $dst" | tee -a "$log_file"
  mkdir -p "$dst"
  echo "Changing to source directory: $src" | tee -a "$log_file"
  cd "$src" || { echo "Error: Source directory not found: $src"; exit 1; }
  echo "Copying files excluding 'docs', 'pkgdown', and all hidden files/directories" | tee -a "$log_file"
  rsync -av --exclude='docs/' --exclude='pkgdown/' --exclude='.*' --exclude='*/.*' --exclude='README.Rmd' --exclude='LICENSE.md' ./ "$dst/"
  echo "Copy Operation Completed Successfully: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$log_file"
} >> "$log_file" 2>&1

echo "Files copied successfully from $src to $dst. Check the log file at $log_file for details."

module load ceuadmin/R
module load mono/5.0.1.1
module load texlive

cd ~/R
export version=$(awk '/Version/{print $2}' ~/pQTLtools/DESCRIPTION)
R CMD build --resave-data --compact-vignettes=both pQTLtools
R CMD INSTALL pQTLtools_${version}.tar.gz
R CMD check --as-cran pQTLtools_${version}.tar.gz
