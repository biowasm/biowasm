#!/bin/bash

DRIVER_CHROME=chromedriver_mac64.zip
DRIVER_FIREFOX=macos

if [[ ! -e chromedriver ]]; then
    version=$(curl http://chromedriver.storage.googleapis.com/LATEST_RELEASE)
    curl -O "http://chromedriver.storage.googleapis.com/${version}/${DRIVER_CHROME}"
    unzip ${DRIVER_CHROME}
    rm ${DRIVER_CHROME}
fi

if [[ ! -e geckodriver ]]; then
    url=$(curl -s https://api.github.com/repos/mozilla/geckodriver/releases/latest | jq -r '.assets[].browser_download_url | select(contains("'$DRIVER_FIREFOX'"))')
    curl -O -L $url
    tar -xvzf $(basename $url)
    rm $(basename $url)
fi
