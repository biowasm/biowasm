// Cloudflare Durable Object that is added to the build via `build.command` in "wrangler.toml".
// Durable Object ID = CDN path, e.g. "seqtk/1.3/seqtk.js"
// Durable Object Names =
// - total: total downloads for one version of a tool
// - YYYY-MM-DD: daily downloads for one version of a tool
export class CDNStats {
	constructor(state) {
		this.state = state;
	}

	// Handle HTTP requests from Cloudflare Worker (called from `[file].js`)
	async fetch(request) {
		const url = new URL(request.url);

		// List all stats, format: { "total": 123, "YYYY-MM-DD": 12 }
		if(url.pathname === "/")
			return new Response(JSON.stringify(
				Object.fromEntries(await this.state.storage.list())
			));

		// Increment daily and total download counts
		if(url.pathname === "/increment") {
			const date = new Date().toISOString().split("T")[0];
			let total = (await this.state.storage.get("total")) || 0;
			let daily = (await this.state.storage.get(date)) || 0;

			// No need for transactions since durable objects feature "automatic write coalescing", so either all writes fail or all succeed
			await this.state.storage.put("total", ++total);
			await this.state.storage.put(date, ++daily);
			return new Response(JSON.stringify({
				total,
				[date]: daily
			}));
		}

		return new Response("Durable Object endpoint not found", { status: 404 });
	}
}
