// GET /cdn/v3/:tool/:version/:file/stats
export async function GET({ request, platform, params }) {
	const path = `${params.tool}/${params.version}/${params.file}`;
	const id = platform.env.stats.idFromName(path);
	const obj = platform.env.stats.get(id);

	// Send HTTP request to Durable Object (needs full path)
	const hostname = new URL(request.url).origin;
	const stats = await (await obj.fetch(hostname)).json();

	return { body: {
		stats
	} };
}
