import { dev } from "$app/environment";
import CONFIG from "@/biowasm.json";

// Download stats, as computed by the cron (see lib/cron.js). Format: {
//   tool: {
//     version: {
//       "YYYY-MM-DD": 123,
//       [...]
//       "total": 123
//     }
//   }
// }
// Server-only: nothing outside this site needs them
export async function getStats(platform) {
	// `platform` is defined even in `vite dev` (adapter-cloudflare emulates it), but with no KV
	// bindings, so the local path has to key off `dev` rather than off `platform`
	if(dev)
		return getMockStats();
	return (await platform.env.CDN.get("STATS", { type: "json" })) || {};
}

// Generate mock stats for local development: a year of daily counts per version
function getMockStats() {
	const stats = {};
	const today = Date.now();
	for(let tool of CONFIG.tools) {
		stats[tool.name] = {};
		for(let version of tool.versions) {
			const counts = { total: 0 };
			for(let day = 365; day > 0; day--) {
				const count = Math.round(100 * Math.random());
				counts[new Date(today - day * 24 * 60 * 60 * 1000).toISOString().slice(0, 10)] = count;
				counts.total += count;
			}
			stats[tool.name][version.version] = counts;
		}
	}

	return stats;
}
