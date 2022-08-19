// Cloudflare Durable Object (added to the Cloudflare Worker logic in `bin/postbuild.sh`).

// Durable Object ID = CDN path, e.g. "seqtk/1.3/seqtk.js"
// Durable Object Names =
// - YYYY-MM-DD: daily downloads for one version of a tool
export class CDNStats {
	constructor(state) {
		this.state = state;
	}

	// Handle HTTP requests from Cloudflare Worker (called from `[file].js`)
	async fetch(request) {
		const url = new URL(request.url);

		// List all stats, format: { "YYYY-MM-DD": 12 }
		if(url.pathname === "/")
			return new Response(JSON.stringify(
				Object.fromEntries(await this.state.storage.list())
			));

		// Increment daily download counts
		if(url.pathname === "/increment") {
			const date = new Date().toISOString().split("T")[0];
			const daily = (await this.state.storage.get(date)) || 0;

			await this.state.storage.put(date, daily + 1);
			return new Response("{}");
		}

		return new Response("Durable Object endpoint not found", { status: 404 });
	}
}
