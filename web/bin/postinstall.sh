#!/bin/bash

# ==============================================================================
# Post-install workarounds (runs after "npm install")
# ==============================================================================

# Set Popper JS to "type=module". Otherwise, we get the following error when using SvelteKit + Sveltestrap:
# 	Cannot use import statement outside a module node_modules/@popperjs/core/dist/esm/popper.js:1
# 	import { popperGenerator, detectOverflow } from "./createPopper.js";
# Using fix from https://github.com/bestguy/sveltestrap/pull/356#issuecomment-1101572342
# Fix from https://github.com/sveltejs/kit/issues/4504#issuecomment-1135338008 did not work
echo "Fix Popper JS module issue..."
FILE=node_modules/@popperjs/core/package.json
jq '.type = "module"' $FILE > $FILE.tmp && mv $FILE.tmp $FILE
