import { isValidTool } from "$lib/utils";
import CONFIG from "@/biowasm.json";

// TODO: GET /api/v3/stats
// TODO: GET /api/v3/:tool/stats
// TODO: GET /api/v3/:tool/:version/stats
export async function GET({ request, platform, params }) {
	let [toolName, version, program] = params.path.split("/", 3);
	toolName ||= null, version ||= null, program ||= null;

	// TODO: Get all stats from KV (more efficient than querying all the tools?)
	if(toolName === null) {
		// console.log({...await platform.env.stats.list()});
		return {
			status: 200,
			body: { stats: "TODO" }
		};
	}

	// Get tool info
	const tool = CONFIG.tools.find(t => t.name === toolName);
	if(!tool)
		return error(params);
	// Get matching version(s)
	const versions = tool.versions.filter(v => version === null || v.version === version);
	if(versions.length === 0)
		return error(params);
	// Get matching program(s)
	const programs = (tool.programs || [toolName]).filter(p => program === null || p === program);
	if(programs.length === 0)
		return error(params);

	// Gather stats for each program
	let stats = {};
	for(let program of programs) {
		stats[program] = {};
		for(let versionInfo of versions)
			stats[program][versionInfo.version] = await getStats(toolName, versionInfo.version, program);
	}

	return {
		status: 200,
		body: { stats }
	};

	// Get stats 
	async function getStats(tool, version, program) {
		// Local dev
		if(platform === undefined)
			return { "2022-01-01": 10, "2022-01-02": 20, "total": 30 };

		// Fetch stats
		const path = `${tool}/${version}/${program}.js`;
		const id = platform.env.stats.idFromName(path);
		const obj = platform.env.stats.get(id);

		// Send HTTP request to Durable Object (needs full path)
		const hostname = new URL(request.url).origin;
		const response = await obj.fetch(hostname);
		return await response.json();
	}	
}

// Return error response
function error(params) {
	return {
		status: 404,
		body: { error: "Could not find tool", params }
	};
}
