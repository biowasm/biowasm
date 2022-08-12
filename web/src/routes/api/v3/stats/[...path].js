import CONFIG from "@/biowasm.json";

// GET /api/v3/stats
// GET /api/v3/stats/:tool
// GET /api/v3/stats/:tool/:version
// GET /api/v3/stats/:tool/:version/:program
export async function GET({ request, platform, params }) {
	let [toolName, versionName, programName] = params.path.split("/", 3);
	toolName ||= null, versionName ||= null, programName ||= toolName;

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
		if(platform === undefined)
			return {
				status: 200,
				body: { stats: getMockStats(toolName, versionName, programName) }
			};

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
	let stats = getMockStats();
	if(platform !== undefined)
		stats = await platform.env.CDN.get("STATS", { type: "json" });

	// Subset stats based on URL parameters
	if(toolName !== null)
		stats = { [toolName]: stats[toolName] };
	if(programName !== null)
		stats[toolName] = { [programName]: stats[toolName][programName] };
	if(versionName !== null)
		stats[toolName][programName] = { [versionName]: stats[toolName][programName][versionName] };

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
			[program]: {
				[version]: stats
			}
		}
	}
}

// Generate mock stats for local development
function getMockStats(tool="samtools", version="1.10", program="samtools") {
	return formatStats({
		"2022-01-01": 10,
		"2022-01-02": 20,
		"total": 30
	}, tool, version, program);
}
