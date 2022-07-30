// Cloudflare Durable Object that is added to the build via `build.command` in "wrangler.toml"
export class CDNStats {
	constructor(state, env) {
		this.state = state;
	}

	// Handle HTTP requests from clients.
	async fetch(request) {
		let url = new URL(request.url);

		// Durable Object storage is automatically cached in-memory, so reading the
		// same key every request is fast. (That said, you could also store the
		// value in a class member if you prefer.)
		let value = (await this.state.storage.get("value")) || 0;

		switch (url.pathname) {
		case "/increment":
			++value;
			break;
		case "/decrement":
			--value;
			break;
		case "/":
			// Just serve the current value.
			break;
		default:
			return new Response("Not found", { status: 404 });
		}

		// You do not have to worry about a concurrent request having modified the
		// value in storage because "input gates" will automatically protect against
		// unwanted concurrency. So, read-modify-write is safe. For more details,
		// refer to: https://blog.cloudflare.com/durable-objects-easy-fast-correct-choose-three/
		await this.state.storage.put("value", value);

		return new Response(value);
	}
}
