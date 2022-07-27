#!/bin/bash

# # Need at least Node v16
# # FIXME: move to deploy.yml
# curl -fsSL https://deb.nodesource.com/setup_16.x | sudo -E bash -
# sudo apt install -y nodejs

# 
npm install
npm run build
cp dist/aioli.js ../build/
