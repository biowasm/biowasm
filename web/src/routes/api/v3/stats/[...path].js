import CONFIG from "@/biowasm.json";
import { getMockStats } from "$lib/utils";

// GET /api/v3/stats
// GET /api/v3/stats/:tool
// GET /api/v3/stats/:tool/:version
// GET /api/v3/stats/:tool/:version/:program
// Format: {
//   stats: {
//     tool: {
//       version: {
//         program: {
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

	// Only get stats from Durable Object if we have a specific tool/version/program
	if(toolName !== null && versionName !== null && programName !== null) {
		// Local dev
		if(platform === undefined) {
			const stats = { "2022-01-01": 10, "2022-01-02": 20, "total": 30 };
			return {
				status: 200,
				body: {
					stats: formatStats(stats, toolName, versionName, programName)
				}
			};
		}

		// Fetch stats
		const path = `${toolName}/${versionName}/${programName}.js`;
		const id = platform.env.stats.idFromName(path);
		const obj = platform.env.stats.get(id);

		// Send HTTP request to Durable Object (`obj.fetch()` needs full path)
		const hostname = new URL(request.url).origin;
		const stats = await (await obj.fetch(hostname)).json();
		return {
			status: 200,
			body: { stats: formatStats(stats, toolName, versionName, programName) }
		};
	}

	// Get stats from KV (faster than querying all Durable Objects)
	const stats = platform === undefined ? getMockStats() : await platform.env.CDN.get("STATS", { type: "json" });

	// Subset stats based on URL parameters
	if(toolName !== null)
		stats = { [toolName]: stats[toolName] };
	if(versionName !== null)
		stats[toolName] = { [versionName]: stats[toolName][versionName] };
	if(programName !== null)
		stats[toolName][versionName] = { [programName]: stats[toolName][versionName][programName] };

	return {
		status: 200,
		body: { stats }
	};
}

// Return error response
function error(params) {
	return {
		status: 404,
		body: { error: "Could not find tool", params }
	};
}

// Format stats object
function formatStats(stats={}, tool="samtools", version="1.10", program="samtools") {
	return {
		[tool]: {
			[version]: {
				[program]: stats
			}
		}
	}
}

