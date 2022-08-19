<script context="module">
import { browser } from "$app/env";
</script>

<script>
import StatsChart from "$components/StatsChart.svelte";
export let data = { series: [] };
</script>

<svelte:head>
	<title>Stats</title>
</svelte:head>

<h4>Stats</h4>

<!-- Show plot (use `browser` check to skip SSR) -->
{#if browser}
	<StatsChart series={data.series} />
{/if}
