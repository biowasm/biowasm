export async function load({ fetch }) {
	// Get stats
	const stats = (await fetch("/api/v3/stats").then(d => d.json())).stats;

	// Aggregate stats across versions
	const statsAgg = {};
	for(let toolName in stats) {
		for(let versionName in stats[toolName]) {
			if(!(toolName in statsAgg))
				statsAgg[toolName] = {};
			for(let date in stats[toolName][versionName]) {
				if(!(date in statsAgg[toolName]))
					statsAgg[toolName][date] = 0;
				statsAgg[toolName][date] += stats[toolName][versionName][date];
			}
		}
	}

	// Convert to Plotly format
	const series = [];
	for(let tool in statsAgg) {
		const x = [];
		const y = [];
		const dates = Object.keys(statsAgg[tool]).filter(d => d !== "total").sort();
		for(let date of dates) {
			x.push(date);
			y.push(statsAgg[tool][date]);
		}
		series.push({ name: tool, x, y });
	}

	return { 
		series
	};
}
