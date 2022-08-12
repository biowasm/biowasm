// Cron that gathers stats from all Durable Objects and stores a summary in a KV
// for fast access. Note that BIOWASM_CONFIG is defined in `bin/postbuild.sh`.
entry_default.scheduled = async (event, env, ctx) => {
	console.log("Starting cron on env =", env.ENVIRONMENT);
	ctx.waitUntil(cron(env));
};

// Cron logic
async function cron(env) {
	const stats = {};

	// Loop through each tool
	for(let tool of BIOWASM_CONFIG.tools) {
		stats[tool.name] = {};

		// Loop through each program
		for(let program of tool.programs) {
			stats[tool.name][program] = {};

			// Loop through each version
			for(let version of tool.versions) {
				// Get Durable Object
				const path = `${tool.name}/${version.version}/${program}.js`;
				const id = env.stats.idFromName(path);
				const obj = env.stats.get(id);

				// Send HTTP request to Durable Object
				const data = await (await obj.fetch(`https://biowasm-v3-${env.ENVIRONMENT}.robert.workers.dev`)).json();
				stats[tool.name][program][version.version] = data;
			}
		}
	}

	// Save in KV store
	await env.CDN.put("STATS", JSON.stringify(stats));
}
