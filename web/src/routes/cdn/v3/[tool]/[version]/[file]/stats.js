import { isValidTool } from "$lib/utils";

// GET /cdn/v3/:tool/:version/:file/stats
export async function GET({ request, platform, params }) {
	// Stop if unrecognized file
	if(!isValidTool(params)) {
		return {
			status: 404,
			body: { error: `Could not find tool`, params }
		};
	}

	const id = platform.env.stats.idFromName(path);
	const obj = platform.env.stats.get(id);

	// Send HTTP request to Durable Object (needs full path)
	const hostname = new URL(request.url).origin;
	const stats = await (await obj.fetch(hostname)).json();

	return { body: {
		stats
	} };
}
