import ASSET_MANIFEST from "../../../biowasm.manifest.json";

export function isValidTool(params) {
	const path = `${params.tool}/${params.version}/${params.file}`;
	return ASSET_MANIFEST[path];
}
