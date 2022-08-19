// On Cloudflare, `ENVIRONMENT` will be stg or prd (see wrangler.toml); locally, it will be undefined
const ENV = typeof ENVIRONMENT === "undefined" ? "stg" : ENVIRONMENT;
const ASSET_MANIFESTS = import.meta.glob("@/biowasm.manifest*.json", { eager: true });
export const ASSET_MANIFEST = ASSET_MANIFESTS[ENV === "prd" ? "../biowasm.manifest.json" : "../biowasm.manifest.stg.json"].default;

// Check if tool is valid
export function isValidTool(params) {
	const path = `${params.tool}/${params.version}/${params.file}`;
	return ASSET_MANIFEST[path];
}

// Generate mock stats for local development
export function getMockStats() {
	const stats = { "2022-07-01": 10, "2022-07-02": 20, "total": 30 };
	return {
		samtools: {
			"1.10": stats
		},
		seqtk: {
			"1.2": stats,
			"1.3": stats
		},
		coreutils: {
			"8.32": stats
		}
	};
}
