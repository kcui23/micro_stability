#!/usr/bin/env bash
set -e

# Set the download URL and file path for the .pkg file
PKG_URL="https://cloud.r-project.org/bin/macosx/big-sur-arm64/base/R-4.4.2-arm64.pkg"
PKG_FILE="r-mac/latest_r.pkg"

mkdir -p r-mac

# Check if the .pkg file already exists
if [ -f "$PKG_FILE" ]; then
  echo "Package file already exists, skipping download."
else
  echo "Downloading package..."
  curl -o "$PKG_FILE" "$PKG_URL"
fi

# Extract the package
cd r-mac
xar -xf latest_r.pkg
rm -r Resources tcltk.pkg texinfo.pkg Distribution latest_r.pkg

# Extract contents of the R framework package
cat R-fw.pkg/Payload | gunzip -dc | cpio -i
mv R.framework/Versions/Current/Resources/* .
rm -r R-fw.pkg R.framework

# Patch the main R script
sed -i.bak '/^R_HOME_DIR=/d' bin/R
sed -i.bak 's;/Library/Frameworks/R.framework/Resources;${R_HOME};g' bin/R
chmod +x bin/R
rm -f bin/R.bak

# Remove unnecessary files
rm -r doc tests
rm -r lib/*.dSYM
