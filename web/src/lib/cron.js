// Cron that gathers stats from all Durable Objects and stores a summary in a KV
// for fast access. Note that BIOWASM_CONFIG is defined in `bin/postbuild.sh`.
entry_default.scheduled = async (event, env, ctx) => {
	ctx.waitUntil(cron(env));
};

// Cron logic
async function cron(env) {
	console.log("Starting cron on env =", env.ENVIRONMENT);

	// Loop through each tool
	const stats = {};
	const latestVersion = {};  // { "seqtk": "1.3", "samtools": "1.10" }
	for(let tool of BIOWASM_CONFIG.tools) {
		stats[tool.name] = {};
		latestVersion[tool.name] = tool.versions.at(-1).version;

		// Loop through each version
		for(let version of tool.versions) {
			stats[tool.name][version.version] = {};

			// Loop through each program
			for(let program of tool.programs) {
				// Get Durable Object
				const path = `${tool.name}/${version.version}/${program}.js`;
				const id = env.stats.idFromName(path);
				const obj = env.stats.get(id);

				// Send HTTP request to Durable Object
				const statsPrgm = await (await obj.fetch(`https://biowasm-v3-${env.ENVIRONMENT}.robert.workers.dev`)).json();
				for(let date in statsPrgm) {
					if(!(date in stats[tool.name][version.version]))
						stats[tool.name][version.version][date] = 0;
					stats[tool.name][version.version][date] += statsPrgm[date];
				}
			}
		}
	}
	// Special case for Aioli: assign CDN v2 stats to 2.x.x versions
	latestVersion["aioli"] = "2.x.x";

	// Bring in stats from v2 CDN
	const statsOld = await env.CDN_V2.get("summary", { type: "json" });
	for(let date in statsOld) {
		for(let toolName in statsOld[date]) {
			// CDN v2 didn't store stats per version so we use latest version
			const versionName = latestVersion[toolName];
			if(!(toolName in stats))
				stats[toolName] = {};
			if(!(versionName in stats[toolName]))
				stats[toolName][versionName] = {};
			if(!(date in stats[toolName][versionName]))
				stats[toolName][versionName][date] = 0
			stats[toolName][versionName][date] += statsOld[date][toolName];
		}
	}

	// Set totals
	for(let toolName in stats) {
		for(let versionName in stats[toolName]) {
			let total = 0;
			for(let date in stats[toolName][versionName])
				if(date !== "total")
					total += stats[toolName][versionName][date];
			stats[toolName][versionName]["total"] = total;
		}
	}

	// Save in KV store
	await env.CDN.put("STATS", JSON.stringify(stats));

	console.log("Done cron on env =", env.ENVIRONMENT);
}
