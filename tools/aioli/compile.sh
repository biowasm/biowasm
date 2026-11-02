#!/bin/bash

# Install a modern Node just for building the Aioli JS library.
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
\. "$HOME/.nvm/nvm.sh"
nvm install 22
nvm use 22

# Build Aioli
npm install
npm run build
cp dist/aioli.js ../build/
