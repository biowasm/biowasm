#!/bin/bash

# Build Aioli
npm install
npm run build
cp dist/aioli.js ../build/
