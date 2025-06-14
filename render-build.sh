#!/usr/bin/env bash
# exit on error
set -o errexit

# Install graphviz
apt-get update && apt-get install -y graphviz

# Install Python dependencies
pip install -r requirements.txt 