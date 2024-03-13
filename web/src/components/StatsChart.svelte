<script>
export let series = [];  // [{ x: [1,2,3], y: [4,5,6] }]
let elPlot;
</script>

<!-- svelte-ignore missing-declaration -->
<svelte:head>
	<script async src="https://cdn.plot.ly/plotly-2.14.0.min.js" on:load={() => {
		// Zoom on data from 3 months ago
		const dateToday = new Date();
		const dateSixMonthsAgo = new Date();
		dateSixMonthsAgo.setMonth(dateToday.getMonth() - 3);

		// Determine initial y-range max
		const yRange = Math.max(...series.map(toolStats => {
			// Find first date in the 6 months time range
			const index = toolStats.x.findIndex(date => new Date(date) > dateSixMonthsAgo)
			const y = toolStats.y.slice(index)

			// Get the max y value for that tool in that date range
			return Math.max(...y);
		}))

		// Plot
		Plotly.newPlot(elPlot, series, {
			title: { text: "Package downloads over time" },
			xaxis: { range: [dateSixMonthsAgo, dateToday].map(d => d.toISOString().slice(0, 10)) },
			yaxis: { range: [0, yRange] },
		});
	}}></script>
</svelte:head>

<div bind:this={elPlot} style="width:100%; height:70vh"></div>
