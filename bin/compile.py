#!/usr/bin/python3.8

import json
import argparse

# State
CONFIG = {}
DIR_BUILD = "build"


def list():
	"""
	List all biowasm packages
	"""
	for tool in CONFIG["tools"]:
		name = tool["name"]
		versions = ", ".join([ v["version"] for v in tool["versions"] ])
		print(f"{name.ljust(10)}\t{versions}")


def compile(tool, versions=[]):
	"""
	Compile one tool, multiple versions
	"""
	# print(f"# Compiling '{tool}' versions '{versions}'")
	tool_info = next(t for t in CONFIG["tools"] if t["name"] == tool)
	files = tool_info["files"] if "files" in tool_info else ["js", "wasm"]
	programs = tool_info["programs"] if "programs" in tool_info else [tool]
	versions = [v for v in tool_info["versions"] if v["version"] in versions] if versions else tool_info["versions"]

	print(f"TOOL={tool} make init")

	for version_info in versions:
		version = version_info["version"]
		branch = version_info["branch"]
		dependencies = version_info["dependencies"] if "dependencies" in version_info else []
		for dependency_info in dependencies:
			compile(dependency_info["name"], [dependency_info["version"]])

		dir_build = f"{DIR_BUILD}/{tool}/{version}"
		print(f"mkdir -p {dir_build}")
		print(f"bin/compile.sh tools/{tool} {branch}")
		for program in programs:
			for file in files:
				print(f"cp tools/{tool}/build/{program}.{file} {dir_build}/")


if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	with open("biowasm.json") as f:
		CONFIG = json.load(f)

	parser = argparse.ArgumentParser()
	parser.add_argument("--list", default=False, action="store_true", help="List tools available")
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
