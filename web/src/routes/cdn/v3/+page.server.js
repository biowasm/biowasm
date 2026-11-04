import { getStats } from "$lib/server/stats";

// Load download stats on the server so the table can be sorted by downloads on first render.
// Only per-tool totals and weekly counts are returned: the full stats (daily counts) are over 1MB
const NB_WEEKS = 52;
const MS_PER_DAY = 24 * 60 * 60 * 1000;

export async function load({ platform }) {
	const stats = await getStats(platform);

	// Map each day of the last 52 complete weeks to its week, so the loop below is a lookup
	// rather than parsing ~80k dates. The current week is left out since its partial count
	// would look like a drop
	const today = new Date();
	const monday = Date.UTC(today.getUTCFullYear(), today.getUTCMonth(), today.getUTCDate()) - ((today.getUTCDay() + 6) % 7) * MS_PER_DAY;
	const dayToWeek = {};
	for(let day = 0; day < NB_WEEKS * 7; day++)
		dayToWeek[new Date(monday - (NB_WEEKS * 7 - day) * MS_PER_DAY).toISOString().slice(0, 10)] = Math.floor(day / 7);

	const downloads = {};
	const weekly = {};
	for(const tool in stats) {
		downloads[tool] = 0;
		weekly[tool] = new Array(NB_WEEKS).fill(0);
		for(const version in stats[tool]) {
			downloads[tool] += stats[tool][version]?.total || 0;

			// Sum daily counts into weeks
			for(const date in stats[tool][version]) {
				const week = dayToWeek[date];
				if(week !== undefined)
					weekly[tool][week] += stats[tool][version][date];
			}
		}
	}

	return { downloads, weekly };
}
