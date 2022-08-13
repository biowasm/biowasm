<script>
export let series = [];  // [{ x: [1,2,3], y: [4,5,6] }]
let elPlot;
</script>

<svelte:head>
	<script async src="https://cdn.plot.ly/plotly-2.14.0.min.js" on:load={() => {
		// Zoom on data from 3 months ago
		let dateToday = new Date();
		let dateSixMonthsAgo = new Date();
		dateSixMonthsAgo.setMonth(dateToday.getMonth() - 3);

		// Plot
		const range = [dateSixMonthsAgo, dateToday].map(d => d.toISOString().slice(0, 10));
		Plotly.newPlot(elPlot, series, { xaxis: { range } });
	}}></script>
</svelte:head>

<div bind:this={elPlot} style="width:100%; height:400px"></div>
