import { json } from "@sveltejs/kit";
import { getAssetFromKV } from "@cloudflare/kv-asset-handler";
import CONFIG from "@/biowasm.json";
const ASSET_MANIFESTS = import.meta.glob("@/biowasm.manifest*.json", { eager: true });

// Settings
const CACHE_CONFIG = {
	browserTTL: 604800,  // 1 week (default: null)
	edgeTTL: 172800,     // 2 days (default: 2 days)
	bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
};
const latestAioliVersion = CONFIG.tools.find(t => t.name === "aioli").versions.at(-1).version;

// GET /cdn/v3/:tool/:version/:file
export async function GET({ request, platform, params }) {
	// Load the right asset manifest file depending on the environment we're in
	const ENV = platform === undefined ? "dev" : platform.env.ENVIRONMENT;
	const ASSET_MANIFEST = ASSET_MANIFESTS[ENV === "prd" ? "../biowasm.manifest.json" : "../biowasm.manifest.stg.json"].default;

	// Stop if unrecognized file
	if(params.version !== "latest" && !(`${params.tool}/${params.version}/${params.file}` in ASSET_MANIFEST)) {
		return json({
			params,
			error: `Could not find tool`
		}, {
			status: 404
		});
	}

	// Local dev
	if(platform === undefined) {
		return json({
			"download": params
		});
	}

	// Artificial /latest route for Aioli
	if(params.tool === "aioli" && params.version === "latest")
		params.version = latestAioliVersion;

	// Download file from Cloudflare Workers
	let response = await getAssetFromKV({
		request,
		// The package `kv-asset-handler` will update Cloudflare Cache after the request is done
		waitUntil: promise => platform.context.waitUntil(promise)
	},
	{
		ASSET_MANIFEST,
		ASSET_NAMESPACE: platform.env.CDN,
		cacheControl: CACHE_CONFIG,
		mapRequestToAsset: request => {
			// Build up URL (can't rely on request.url containing the tool/version/file since it can be called from /latest/)
			const url = `${new URL(request.url).origin}/${params.tool}/${params.version}/${params.file}`;
			return new Request(url, request);
		},
	});

	// Log event only after return file to user
	const path = `${params.tool}/${params.version}/${params.file}`;
	platform.context.waitUntil(logEvent({ request, platform, path }));

	// Enable CORS
	response.headers.set("Access-Control-Allow-Origin", "*");
	return response;
}

async function logEvent({ request, platform, path }) {
	const hostname = new URL(request.url).origin;

	// Get (or create) a durable object based on a valid path
	const id = platform.env.stats.idFromName(path);
	const obj = platform.env.stats.get(id);

	// Send HTTP request to Durable Object (needs full URL, just "/increment" doesn't work)
	await obj.fetch(`${hostname}/increment`);
}
