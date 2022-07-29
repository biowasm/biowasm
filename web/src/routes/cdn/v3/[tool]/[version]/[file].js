import { getAssetFromKV, mapRequestToAsset } from "@cloudflare/kv-asset-handler";
import ASSET_MANIFEST from "../../../../../../../biowasm.manifest.json";

// Settings
const URL_CDN = "cdn/v3";
const CACHE_CONFIG = {
	browserTTL: 604800,  // 1 week (default: null)
	edgeTTL: 172800,     // 2 days (default: 2 days)
	bypassCache: false   // Do not bypass Cloudflare's cache (default: false)
};

// GET /cdn/v3/:tool/:version/:file
export async function GET({ request, platform, params }) {
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
			return mapRequestToAsset(new Request(url, request));
		},
	});

	// Enable CORS
	response.headers.set("Access-Control-Allow-Origin", "*");
	return response;
}
