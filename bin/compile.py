#!/usr/bin/python3
# This script compiles tools of interest to WebAssembly and regenerates "biowasm.manifest.json"

import json
import base64
import argparse
import subprocess
from pathlib import Path

# Config
DIR_ROOT = Path(__file__).parent / "../"
DIR_BUILD = (DIR_ROOT / "build").resolve()
DIR_MANIFEST = (DIR_ROOT / "biowasm.manifest.json").resolve()
DIR_MANIFEST_TEMP = (DIR_BUILD / "manifest.tmp").resolve()
DIR_CF_UPLOAD = (DIR_BUILD / "cf_kv_upload.json").resolve()
DIR_CONFIG = (DIR_ROOT / "biowasm.json").resolve()
COLOR_GREEN = "\033[1;33m"
COLOR_OFF = "\033[0m"

# Global state
CONFIG = {}
MANIFEST = {}
LEVEL = 0


def list():
	"""
	List all biowasm packages
	"""
	for tool in CONFIG["tools"]:
		name = tool["name"]
		versions = ", ".join([ v["version"] for v in tool["versions"] ])
		print(f"{name.ljust(10)}\t{versions}")


def exec(cmd):
	"""
	Execute a Bash command or just print it on screen if dry run enabled
	"""
	cmd_dry = f"{COLOR_GREEN}{'    ' * LEVEL}{cmd}{COLOR_OFF}"
	print(cmd_dry)

	if not args.dry_run:
		# Stream output to screen as we get it
		with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True) as p:
			for line in p.stdout:
				print(line, end='')
		if p.returncode != 0:
			print("Return code not 0: ", subprocess.CalledProcessError(p.returncode, p.args))


def get_file_contents(path):
	"""
	Get file contents as text if .js, or base64 if .wasm/.data
	"""
	if args.dry_run:
		return f"<path>"

	if path.endswith(".wasm") or path.endswith(".data"):
		with open(path, "rb") as f:
			return base64.b64encode(f.read()).decode("utf-8"), True
	elif path.endswith(".js"):
		with open(path, "r") as f:
			return f.read(), False
	else:
		raise f"Unexpected file extension for file '{path}'."


def compile(tool, versions=[], level=0):
	"""
	Compile one tool, multiple versions
	"""
	global LEVEL
	LEVEL = level

	# Get tool info
	tool_info = next(t for t in CONFIG["tools"] if t["name"] == tool)
	files = tool_info["files"] if "files" in tool_info else ["js", "wasm"]
	programs = tool_info["programs"] if "programs" in tool_info else [tool]
	versions = [v for v in tool_info["versions"] if v["version"] in versions] if versions else tool_info["versions"]
	tool_git_path = f"tools/{tool}/src/"

	# Validate
	if len(versions) == 0:
		print(f"No valid versions found in biowasm.json for {tool}. Make sure versions don't start with 'v'.")
		exit(1)

	# Init repo
	exec(f"git submodule update --init --recursive {tool_git_path} && git submodule status {tool_git_path}")

	# Compile each version and its dependencies
	for version_info in versions:
		version = version_info["version"]
		branch = version_info["branch"]
		dependencies = version_info["dependencies"] if "dependencies" in version_info else []
		for dependency_info in dependencies:
			compile(tool=dependency_info["name"], versions=[dependency_info["version"]], level=level + 1)
		LEVEL = level

		dir_build = f"{DIR_BUILD}/{tool}/{version}"
		exec(f"mkdir -p {dir_build}")
		exec(f"bin/compile.sh tools/{tool} {branch}")
		for program in programs:
			for file in files:
				exec(f"cp tools/{tool}/build/{program}.{file} {dir_build}/")
				exec(f"md5sum {dir_build}/{program}.{file} | sed 's|{DIR_BUILD}/||' >> {DIR_MANIFEST_TEMP}")


def generate_manifests():
	"""
	Regenerate manifest files
	"""
	if args.dry_run:
		print("<update manifest>")
		return

	# Update the MANIFEST where needed
	to_upload = []
	with open(DIR_MANIFEST_TEMP) as f:
		for row in f.readlines():
			hash, path = row.rstrip().split("  ")

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
	with open(DIR_CF_UPLOAD, "w", encoding="utf-8") as f:
		json.dump(to_upload, f, ensure_ascii=False)
	with open(DIR_MANIFEST, "w", encoding="utf-8") as f:
		json.dump(MANIFEST, f, ensure_ascii=False, sort_keys=True, indent=2)
	exec(f"rm {DIR_MANIFEST_TEMP}")


if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--list", default=False, action="store_true", help="List tools available")
	parser.add_argument("--dry-run", default=False, action="store_true", help="Dry run")
	parser.add_argument("--tools", type=str, help="Tool to compile to WebAssembly")
	parser.add_argument("--versions", type=str, help="Tool version(s)")
	parser.add_argument("--env", type=str, help="Environment (stg, prd)", default="stg")
	args = parser.parse_args()

	# Load configs
	if args.env == "stg":
		DIR_MANIFEST = Path(str(DIR_MANIFEST).replace('.json', '.stg.json'))
		DIR_MANIFEST_TEMP = Path(str(DIR_MANIFEST_TEMP).replace('.tmp', '.stg.tmp'))

	if DIR_MANIFEST_TEMP.is_file():
		print("Manifest file already exists from previous run. Delete it first.")
		exit(1)

	with open(DIR_CONFIG) as f:
		CONFIG = json.load(f)
	with open(DIR_MANIFEST) as f:
		MANIFEST = json.load(f)

	# Process CLI parameters
	if args.list:
		list()
	elif args.tools:
		# Compile all tools
		tools = args.tools.split(",") if args.tools != "all" else [ t["name"] for t in CONFIG["tools"] ]
		for tool_name in tools:
			compile(tool_name, args.versions.split(",") if args.versions is not None else [])
		# Regenerate manifest files
		generate_manifests()
	else:
		parser.print_usage()
