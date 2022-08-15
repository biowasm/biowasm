<script>
export let series = [];  // [{ x: [1,2,3], y: [4,5,6] }]
let elPlot;
</script>

<!-- svelte-ignore missing-declaration -->
<svelte:head>
	<script async src="https://cdn.plot.ly/plotly-2.14.0.min.js" on:load={() => {
		// Zoom on data from 3 months ago
		let dateToday = new Date();
		let dateSixMonthsAgo = new Date();
		dateSixMonthsAgo.setMonth(dateToday.getMonth() - 3);

		// Plot
		Plotly.newPlot(elPlot, series, {
			title: { text: "Package downloads over time" },
			xaxis: { range: [dateSixMonthsAgo, dateToday].map(d => d.toISOString().slice(0, 10)) }
		});
	}}></script>
</svelte:head>

<div bind:this={elPlot} style="width:100%; height:70vh"></div>
