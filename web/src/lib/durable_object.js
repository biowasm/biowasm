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
		const date = new Date().toISOString().split("T")[0];
		let total = (await this.state.storage.get("total")) || 0;
		let daily = (await this.state.storage.get(date)) || 0;

		const url = new URL(request.url);
		switch (url.pathname) {
			// Increment daily and total download counts
			case "/increment":
				total++;
				daily++;
				await this.state.storage.put("total", total);
				await this.state.storage.put(date, daily);
				break;

			// Fetch values
			case "/":
				break;

			default:
				return new Response("Durable Object endpoint not found", { status: 404 });
		}

		return new Response(JSON.stringify({
			total,
			daily
		}));
	}
}
