import { getAssetFromKV } from "@cloudflare/kv-asset-handler";
import CONFIG from "../../../../../../../biowasm.json";
import ASSET_MANIFEST from "../../../../../../../biowasm.manifest.json";

// Settings
const URL_CDN = "/cdn/v3";
const CACHE_CONFIG = {
	browserTTL: 604800,  // 1 week (default: null)
	edgeTTL: 172800,     // 2 days (default: 2 days)
	bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
};
const latestAioliVersion = CONFIG.tools.find(t => t.name === "aioli").versions.at(-1).version;

// GET /cdn/v3/:tool/:version/:file
export async function GET({ request, platform, params }) {
	// Unrecognized file
	const path = `${params.tool}/${params.version}/${params.file}`;
	if(!ASSET_MANIFEST[path] && params.version !== "latest") {
		return {
			status: 404,
			body: {
				error: `Could not find tool`,
				params,
		} };
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
		waitUntil: promise => platform.context.waitUntil(promise)
	},
	{
		ASSET_MANIFEST,
		ASSET_NAMESPACE: platform.env.CDN,
		cacheControl: CACHE_CONFIG,
		mapRequestToAsset: request => {
			let url = request.url;
			url = url.replace(URL_CDN, "");
			// Artificial /latest route for Aioli
			if(params.tool === "aioli" && params.version === "latest")
				url = url.replace("latest", latestAioliVersion);
			return new Request(url, request);
		},
	});

	// Enable CORS
	response.headers.set("Access-Control-Allow-Origin", "*");
	return response;
}
