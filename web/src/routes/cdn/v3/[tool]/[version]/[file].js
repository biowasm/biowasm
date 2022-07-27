import { getAssetFromKV, mapRequestToAsset } from "@cloudflare/kv-asset-handler";

export async function GET({ params, request, platform }) {
	const ASSET_MANIFEST = {
		"cdn/v3/samtools/1.10/samtools.json": "test123"
	};

	console.log("platform =", platform);

	return getAssetFromKV(
		{
			request,
			waitUntil: promise => platform.context.waitUntil(promise)
		},
		{
			ASSET_MANIFEST,
			ASSET_NAMESPACE: platform.env.CDN,
			mapRequestToAsset: request => {
				let url = request.url
				url = url.replace("/docs", "").replace(/^\/+/, "")
				return mapRequestToAsset(new Request(url, request))
			}
		}
	);
}
