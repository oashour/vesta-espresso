#!/bin/bash

# This script is used to run the application on Mac OS X.
echo "Running Vesta Espresso on Mac OS X..."
echo "Arguments are $@"

"$PWD/vesta-espresso.dist/vesta-espresso" "$@"