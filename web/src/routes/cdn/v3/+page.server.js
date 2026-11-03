// Load download totals on the server so the table can be sorted by downloads on first render.
// Only the per-tool totals are returned: the full stats (daily counts) are over 1MB
export async function load({ fetch }) {
	const stats = (await (await fetch("/api/v3/stats")).json()).stats || {};

	const downloads = {};
	for(const tool in stats) {
		downloads[tool] = 0;
		for(const version in stats[tool])
			downloads[tool] += stats[tool][version]?.total || 0;
	}

	return { downloads };
}
