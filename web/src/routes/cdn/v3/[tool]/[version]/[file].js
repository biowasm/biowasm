import { getAssetFromKV } from "@cloudflare/kv-asset-handler";
import CONFIG from "../../../../../../../biowasm.json";
import ASSET_MANIFEST from "../../../../../../../biowasm.manifest.json";
import { isValidTool } from "$lib/utils";

// Settings
const CACHE_CONFIG = {
	browserTTL: 604800,  // 1 week (default: null)
	edgeTTL: 172800,     // 2 days (default: 2 days)
	bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
};
const latestAioliVersion = CONFIG.tools.find(t => t.name === "aioli").versions.at(-1).version;

// GET /cdn/v3/:tool/:version/:file
export async function GET({ request, platform, params }) {
	// Stop if unrecognized file
	if(params.version !== "latest" && !isValidTool(params)) {
		return {
			status: 404,
			body: { error: `Could not find tool`, params }
		};
	}

	// Local dev
	if(platform === undefined) {
		return { body: {
			"download": params
		} };
	}

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
			let url = request.url;
			url = url.replace(CONFIG.url, "");
			// Artificial /latest route for Aioli
			if(params.tool === "aioli" && params.version === "latest")
				url = url.replace("latest", latestAioliVersion);
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
	const id = platform.env.stats.idFromName(path);
	const obj = platform.env.stats.get(id);

	// Send HTTP request to Durable Object (needs full path)
	const hostname = new URL(request.url).origin;
	const counts = await (await obj.fetch(`${hostname}/increment`)).json();
	console.log("Counts:", counts);
}
