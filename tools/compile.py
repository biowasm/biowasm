#!/usr/bin/python3.8

import json
import argparse

# State
CONFIG = {}


def list():
	"""
	
	"""
	print(CONFIG)


def compile(tools, versions=None):
	"""
	"""
	for tool in tools:
		tool = next(t for t in CONFIG["tools"] if t["name"] == tool)
		print(f"make init {tool['name']}")
		if versions != ['']:
			tool["versions"] = [ v for v in tool["versions"] if v["version"] in versions ]
		for version in tool["versions"]:
			print(f"tools/compile.sh tools/{tool['name']} {version['branch']}")


if __name__ == "__main__":
	"""
	Parse user arguments
	"""
	with open("biowasm.json") as f:
		CONFIG = json.load(f)

	parser = argparse.ArgumentParser()
	parser.add_argument("--list", default=False, action="store_true", help="List tools available")
	parser.add_argument("--tools", type=str, help="Tool to compile to WebAssembly")
	parser.add_argument("--versions", default="", type=str, help="Tool version(s)")
	args = parser.parse_args()

	if args.list:
		list()
	elif args.tools:
		compile(args.tools.split(","), args.versions.split(","))
	else:
		parser.print_usage()
