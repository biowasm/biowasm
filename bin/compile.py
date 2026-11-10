#!/usr/bin/python3
# This script compiles tools of interest to WebAssembly. Run bin/manifest.py afterwards to turn
# build/ into the manifest and the Cloudflare KV payloads.

import json
import argparse
import subprocess
from pathlib import Path

# Config
DIR_ROOT = Path(__file__).parent / "../"
DIR_BUILD = (DIR_ROOT / "build").resolve()
# Hashes accumulate here across runs so tools can be compiled separately and deployed
# together; bin/manifest.py consumes the file and deletes it
DIR_MANIFEST_TEMP = (DIR_BUILD / "manifest.tmp").resolve()
DIR_CONFIG = (DIR_ROOT / "biowasm.json").resolve()
COLOR_GREEN = "\033[1;33m"
COLOR_OFF = "\033[0m"

# Global state
CONFIG = {}
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
		# Stop here: carrying on would build whatever source happens to be checked out
		if p.returncode != 0:
			print("Return code not 0: ", subprocess.CalledProcessError(p.returncode, p.args))
			exit(1)


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
	# No URL means the source is in this repo already, like "base"
	tool_url = tool_info.get("url")
	# Repos cloned next to src/, e.g. HMMER needs Easel and SIMDe
	tool_sources = tool_info.get("sources", [])

	# Validate
	if len(versions) == 0:
		print(f"No valid versions found in biowasm.json for {tool}. Make sure versions don't start with 'v'.")
		exit(1)

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

		# Grab the source at this version's tag
		if tool_url:
			exec(f"bin/fetch.sh tools/{tool}/src {tool_url} {branch}")
		for source in tool_sources:
			exec(f"bin/fetch.sh tools/{tool}/{source['path']} {source['url']} {source['branch']}")

		exec(f"bin/compile.sh tools/{tool} {branch}")
		for program in programs:
			for file in files:
				exec(f"cp tools/{tool}/build/{program}.{file} {dir_build}/")
				exec(f"md5sum {dir_build}/{program}.{file} | sed 's|{DIR_BUILD}/||' >> {DIR_MANIFEST_TEMP}")


if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--list", default=False, action="store_true", help="List tools available")
	parser.add_argument("--dry-run", default=False, action="store_true", help="Dry run")
	parser.add_argument("--tools", type=str, help="Tool to compile to WebAssembly")
	parser.add_argument("--versions", type=str, help="Tool version(s), comma-separated; 'all' for every version of each tool", default="all")
	parser.add_argument("--env", type=str, help="Environment (stg, prd)", default="stg")
	args = parser.parse_args()

	# Load configs
	if args.env == "stg":
		DIR_MANIFEST_TEMP = Path(str(DIR_MANIFEST_TEMP).replace('.tmp', '.stg.tmp'))

	with open(DIR_CONFIG) as f:
		CONFIG = json.load(f)

	# Process CLI parameters
	if args.list:
		list()
	elif args.tools:
		# Compile all tools
		tools = args.tools.split(",") if args.tools != "all" else [ t["name"] for t in CONFIG["tools"] ]
		for tool_name in tools:
			compile(tool_name, args.versions.split(",") if args.versions not in (None, "", "all") else [])
		print(f"Done. Run `bin/manifest.py --env {args.env}` to update the manifest.")
	else:
		parser.print_usage()
