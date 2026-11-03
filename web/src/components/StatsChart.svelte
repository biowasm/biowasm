<script>
import { onMount } from "svelte";

export let series = [];  // [{ x: [1,2,3], y: [4,5,6] }]
let elPlot;

function plot() {
	// Open on the last month of data (the full series is still plotted, so
	// zooming out shows the rest)
	const dateToday = new Date();
	const dateStart = new Date();
	dateStart.setMonth(dateToday.getMonth() - 1);

	// Determine initial y-range max
	const yRange = Math.max(...series.map(toolStats => {
		// Find the first date inside the visible range
		const index = toolStats.x.findIndex(date => new Date(date) > dateStart)
		const y = toolStats.y.slice(index)

		// Get the max y value for that tool in that date range
		return Math.max(...y);
	}))

	// Plot
	window.Plotly.newPlot(elPlot, series, {
		title: { text: "Package downloads over time" },
		xaxis: { range: [dateStart, dateToday].map(d => d.toISOString().slice(0, 10)) },
		yaxis: { range: [0, yRange] },
	});
}

// If Plotly loaded before this component mounted (e.g. from cache), its load
// event already fired, so plot now instead of waiting for it
onMount(() => {
	if(window.Plotly)
		plot();
});
</script>

<svelte:head>
	<script async src="https://cdn.plot.ly/plotly-2.14.0.min.js" on:load={plot}></script>
</svelte:head>

<div bind:this={elPlot} style="width:100%; height:70vh"></div>
