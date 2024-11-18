#!/bin/bash
rm -rf ./dist
rm -rf ./r-mac
rm -rf /Applications/MicroStability.app
./r-setup.sh

npm run build:mac && open ./dist/micro-stability-1.0.0.dmg