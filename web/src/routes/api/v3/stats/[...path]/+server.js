import { json } from "@sveltejs/kit";
import { dev } from "$app/environment";
import CONFIG from "@/biowasm.json";

const STG = "https://stg.biowasm.com";

// GET /api/v3/stats
// GET /api/v3/stats/:tool
// GET /api/v3/stats/:tool/:version
// GET /api/v3/stats/:tool/:version/:program
// Format: {
//   stats: {
//     tool: {
//       version: {
//         program: {  // <-- only if `program` is specified in URL
//           "YY-MM-DD": 123,
//           [...]
//           "total": 123
//         }
//       }
//     }
//   }
// }
export async function GET({ request, platform, params }) {
	let [toolName, versionName, programName] = params.path.split("/", 3);
	toolName ||= null, versionName ||= null, programName ||= null;

	// Input validation
	if(toolName !== null) {
		const tool = CONFIG.tools.find(t => t.name === toolName);
		if(!tool)
			return error(params);

		// Get matching version(s)
		if(versionName !== null) {
			const version = tool.versions.find(v => v.version === versionName);
			if(!version)
				return error(params);

			// Get matching program(s)
			if(programName !== null) {
				const program = tool.programs.find(p => p === programName);
				if(!program)
					return error(params);
			}
		}
	}

	// Get stats from KV (faster than querying all Durable Objects). Note that `platform` is
	// defined even in `vite dev` (adapter-cloudflare emulates it), but with no KV bindings,
	// so the local path has to key off `dev` rather than off `platform`
	let stats = (dev ? await getLocalStats() : await platform.env.CDN.get("STATS", { type: "json" })) || {};

	// Subset stats based on URL parameters. A tool that is in biowasm.json but has no
	// downloads yet isn't in STATS at all, so don't assume it's there
	if(toolName !== null)
		stats = { [toolName]: stats[toolName] || {} };
	if(versionName !== null)
		stats[toolName] = { [versionName]: stats[toolName][versionName] || {} };

	return json({ stats });
}

// Return error response
function error(params) {
	return json({
		params,
		error: "Could not find tool"
	}, {
		status: 404,
	});
}

// Real stats for local development, borrowed from stg so the numbers on the page are
// plausible instead of random or all-zero. Falls back to mock data when stg is unreachable
async function getLocalStats() {
	try {
		const response = await fetch(`${STG}/api/v3/stats`);
		if(!response.ok)
			throw new Error(`stg returned ${response.status}`);
		return (await response.json()).stats;
	} catch (error) {
		console.warn(`Could not reach ${STG}, falling back to mock stats:`, error.message);
		return getMockStats();
	}
}

// Generate mock stats for local development
function getMockStats() {
	const stats = {};
	for(let tool of CONFIG.tools) {
		stats[tool.name] = {};
		for(let version of tool.versions)
			stats[tool.name][version.version] = {
				"2022-07-01": Math.round(1000 * Math.random()),
				"2022-07-02": Math.round(2000 * Math.random()),
				"total": Math.round(3000 * Math.random())
			};
	}

	return stats;
}
