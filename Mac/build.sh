#!/bin/bash

PYTHON="/usr/local/bin/python3.11"
TOPDIR="$(dirname "$(pwd)")"
BUILD_DIR="build/v$1"

echo "Building VESTA (ESPRESSO) v$1 for Mac..."
echo "Build directory: $BUILD_DIR"
echo "Removing old build..."
rm -rf "$BUILD_DIR"
echo "Making directory structure..."
mkdir -p "$BUILD_DIR/dmg"
mkdir -p "$BUILD_DIR/nuitka"

echo "Making venv..."
cd "$BUILD_DIR/nuitka"
$PYTHON -m venv "venv/"
source "venv/bin/activate"
pip install $TOPDIR nuitka
# This symlink just helps keep the file names clean, 
# not sure how to adjust the nuitka command to do the names right.
ln -sf "$TOPDIR/vesta_espresso/cli.py" "$TOPDIR/vesta_espresso/vesta-espresso.py"
python -m nuitka --standalone "$TOPDIR/vesta_espresso/vesta-espresso.py" -o vesta-espresso
rm "$TOPDIR/vesta_espresso/vesta-espresso.py"
cd -

echo "Calling Platypus..."
platypus \
  --load-profile "vesta-espresso.platypus" \
  --app-version "$1" \
  --bundled-file "$BUILD_DIR/nuitka/vesta-espresso.dist" \
  "$BUILD_DIR/dmg/VESTA (ESPRESSO).app"

echo "Calling create-dmg..."
$HOME/code/create-dmg/create-dmg \
  --volname "VESTA (ESPRESSO)" \
  --volicon "vesta-espresso.icns" \
  --window-pos 200 120 \
  --window-size 650 300 \
  --icon-size 144 \
  --icon "VESTA (ESPRESSO).app" 110 100 \
  --hide-extension "VESTA (ESPRESSO).app" \
  --background "dmg_background.png" \
  --app-drop-link 495 100 \
  --install-destination "/Applications/VESTA" \
  "build/v$1/vesta-espresso-v$1-mac.dmg" \
  "build/v$1/dmg"

echo "Deleting temporary files..."
rm -rf "$BUILD_DIR/nuitka"
rm -rf "$BUILD_DIR/dmg"
