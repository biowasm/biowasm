#!/usr/bin/python3
# This script regenerates "biowasm.manifest.json" and the Cloudflare KV payloads from the
# build/ folder that compile.py generated. Kept separate from compile.py because the manifest
# records what is live in KV, which is only true once the deploy succeeds.

import os
import json
import base64
import argparse
from pathlib import Path

# Config
DIR_ROOT = Path(__file__).parent / "../"
DIR_BUILD = (DIR_ROOT / "build").resolve()
DIR_MANIFEST = (DIR_ROOT / "biowasm.manifest.json").resolve()
# Byte sizes for the website to show. Kept out of the manifest because that one is handed
# straight to getAssetFromKV, which expects a flat path -> key map and nothing else. Not split
# by environment, unlike the manifest: a file's size comes from the build, and a per-env file
# would start out empty and blank the size of every tool that run didn't recompile
DIR_SIZES = (DIR_ROOT / "biowasm.sizes.json").resolve()
DIR_MANIFEST_TEMP = (DIR_BUILD / "manifest.tmp").resolve()
DIR_CF_UPLOAD = (DIR_BUILD / "cf_kv_upload").resolve()  # suffixed with ".<n>.json"
# Cloudflare rejects bulk KV requests over 100MB, and `wrangler kv bulk put` only splits
# on key count, so split on size here instead
CF_UPLOAD_MAX_BYTES = 80 * 1000 * 1000

# Global state
MANIFEST = {}
SIZES = {}


def get_file_contents(path):
	"""
	Get file contents as text if .js, or base64 if .wasm/.data
	"""
	if path.endswith(".wasm") or path.endswith(".data"):
		with open(path, "rb") as f:
			return base64.b64encode(f.read()).decode("utf-8"), True
	elif path.endswith(".js"):
		with open(path, "r") as f:
			return f.read(), False
	else:
		raise ValueError(f"Unexpected file extension for file '{path}'.")


def write_cf_upload(keyvalues, n):
	"""
	Write one batch of key/value pairs for `wrangler kv bulk put`
	"""
	path = f"{DIR_CF_UPLOAD}.{n}.json"
	print(f"Writing {len(keyvalues)} KVs to {path}...")
	with open(path, "w", encoding="utf-8") as f:
		json.dump(keyvalues, f, ensure_ascii=False)


def generate_manifests():
	"""
	Regenerate manifest files
	"""
	# A path shows up more than once when a tool was rebuilt, whether that's a failed attempt
	# followed by a good one or two compile.py runs days apart, so the last hash wins
	hashes = {}
	with open(DIR_MANIFEST_TEMP) as f:
		for row in f.readlines():
			hash, path = row.rstrip().split("  ")
			hashes[path] = hash

	# Update the MANIFEST where needed
	to_upload = []
	for path, hash in hashes.items():
		SIZES[path] = os.path.getsize(DIR_BUILD / path)

		# Prepare the KV pair for Cloudflare Workers KV
		kv_key = f"{path}:{hash}"
		kv_value, kv_base64 = get_file_contents(str(DIR_BUILD / path))

		# If adding new file, or hash of a file has changed, update the manifest
		if path not in MANIFEST or MANIFEST[path] != kv_key:
			print(f"Adding {path} to list of KVs to upload...")
			MANIFEST[path] = kv_key
			to_upload.append({
				"key": kv_key,
				"value": kv_value,
				"base64": kv_base64
			})

	# Save manifest JSON files
	for path in DIR_BUILD.glob("cf_kv_upload.*.json"):
		path.unlink()
	chunk, chunk_bytes, n = [], 0, 0
	for kv in to_upload:
		kv_bytes = len(kv["key"]) + len(kv["value"])
		if chunk and chunk_bytes + kv_bytes > CF_UPLOAD_MAX_BYTES:
			n += 1
			write_cf_upload(chunk, n)
			chunk, chunk_bytes = [], 0
		chunk.append(kv)
		chunk_bytes += kv_bytes
	if chunk:
		write_cf_upload(chunk, n + 1)
	with open(DIR_MANIFEST, "w", encoding="utf-8") as f:
		json.dump(MANIFEST, f, ensure_ascii=False, sort_keys=True, indent=2)
	with open(DIR_SIZES, "w", encoding="utf-8") as f:
		json.dump(SIZES, f, ensure_ascii=False, sort_keys=True, indent=2)
	DIR_MANIFEST_TEMP.unlink()


if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--env", type=str, help="Environment (stg, prd)", default="stg")
	args = parser.parse_args()

	if args.env == "stg":
		DIR_MANIFEST = Path(str(DIR_MANIFEST).replace('.json', '.stg.json'))
		DIR_MANIFEST_TEMP = Path(str(DIR_MANIFEST_TEMP).replace('.tmp', '.stg.tmp'))

	# Without it there are no hashes to fold in, which would silently deploy an unchanged manifest
	if not DIR_MANIFEST_TEMP.is_file():
		print(f"{DIR_MANIFEST_TEMP} not found. Run bin/compile.py first.")
		exit(1)

	with open(DIR_MANIFEST) as f:
		MANIFEST = json.load(f)
	# Sizes of tools that weren't recompiled this run have to survive it
	if DIR_SIZES.is_file():
		with open(DIR_SIZES) as f:
			SIZES = json.load(f)

	generate_manifests()
