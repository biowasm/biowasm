#!/usr/bin/python3.8

import os
import json
import argparse
import subprocess
from pathlib import Path

# Config
DIR_BUILD = "build"
DIR_CONFIG = (Path(__file__).parent / "../biowasm.json").resolve()
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


def exec(cmd):
	cmd_dry = f"{COLOR_GREEN}{'    ' * LEVEL}{cmd}{COLOR_OFF}"
	print(cmd_dry)

	if not args.dry_run:
		# Not safe to use shell=True generally but this is only running trusted code
		output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode().strip()
		if output:
			print(output)

if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	with open(DIR_CONFIG) as f:
		CONFIG = json.load(f)

	parser = argparse.ArgumentParser()
	parser.add_argument("--list", default=False, action="store_true", help="List tools available")
	parser.add_argument("--dry-run", default=False, action="store_true", help="Dry run")
	parser.add_argument("--tools", type=str, help="Tool to compile to WebAssembly")
	parser.add_argument("--versions", type=str, help="Tool version(s)")
	args = parser.parse_args()

	if args.list:
		list()
	elif args.tools:
		for tool_name in args.tools.split(","):
			compile(tool_name, args.versions.split(",") if args.versions is not None else [])
	else:
		parser.print_usage()
