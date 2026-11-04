import { getStats } from "$lib/server/stats";

// Download count for the version page's badge
export async function load({ params, platform }) {
	if(!params.tool || !params.version)
		return { nbDownloads: null };

	const stats = await getStats(platform);
	return { nbDownloads: stats[params.tool]?.[params.version]?.total || 0 };
}
