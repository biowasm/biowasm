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
	for(let tool of BIOWASM_CONFIG.tools) {
		stats[tool.name] = {};

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
				const data = await (await obj.fetch(`https://biowasm-v3-${env.ENVIRONMENT}.robert.workers.dev`)).json();
				stats[tool.name][version.version][program] = data;
			}
		}
	}

	// Save in KV store
	await env.CDN.put("STATS", JSON.stringify(stats));

	console.log("Done cron on env =", env.ENVIRONMENT);
}
