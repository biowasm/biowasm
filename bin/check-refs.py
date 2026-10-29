#!/usr/bin/env python3
"""
Check that every branch/tag in biowasm.json actually exists in the tool's repo.
"""
import json, subprocess, sys

import os, re

def exists(url, ref, local=None):
    """Return "ok", "missing", or "sha" (a commit, which ls-remote cannot resolve)."""
    out = subprocess.run(["git", "ls-remote", url, ref, f"refs/tags/{ref}", f"refs/heads/{ref}"],
                         capture_output=True, text=True, timeout=60).stdout
    if out.strip():
        return "ok"
    # A bare commit SHA is a valid checkout target but is never listed by ls-remote
    if re.fullmatch(r"[0-9a-f]{7,40}", ref):
        if local and os.path.isdir(local):
            r = subprocess.run(["git", "-C", local, "cat-file", "-e", f"{ref}^{{commit}}"],
                               capture_output=True)
            return "ok" if r.returncode == 0 else "missing"
        # Can't check an abbreviated SHA remotely: ls-remote only lists refs, and a
        # shallow fetch needs the full 40 characters
        return "sha"
    return "missing"

bad = 0
unverified = []
for tool in json.load(open("biowasm.json"))["tools"]:
    url = tool.get("url")
    if not url:
        continue
    for v in tool["versions"]:
        r = exists(url, v["branch"], f"tools/{tool['name']}/src")
        if r == "missing":
            bad += 1
            print(f"MISSING    {tool['name']} {v['version']}: ref '{v['branch']}' not in {url}")
        elif r == "sha":
            unverified.append(f"{tool['name']} {v['version']} ({v['branch']})")
    for s in tool.get("sources", []):
        if exists(s["url"], s["branch"], f"tools/{tool['name']}/{s['path']}") == "missing":
            bad += 1
            print(f"MISSING    {tool['name']} source {s['path']}: ref '{s['branch']}' not in {s['url']}")

if unverified:
    print(f"\n{len(unverified)} commit SHA(s), valid but only checkable once cloned:")
    for u in unverified: print(f"  {u}")
print(f"\n{bad} bad ref(s)" if bad else "\nAll named refs resolve.")
sys.exit(1 if bad else 0)
