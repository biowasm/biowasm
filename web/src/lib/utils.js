import ASSET_MANIFEST from "@/biowasm.manifest.json";

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
			"1.10": {
				samtools: stats,
			}
		},
		seqtk: {
			"1.2": {
				seqtk: stats
			},
			"1.3": {
				seqtk: stats
			}
		},
		coreutils: {
			"8.32": {
				head: stats,
				tail: stats,
				wc: stats
			}
		}
	};
}
