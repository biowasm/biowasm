<script>
// Tiny line chart of weekly downloads, scaled to its own max
export let values = [];
export let width = 100;
export let height = 24;

$: max = Math.max(1, ...values);
$: points = values.map((v, i) => [i / Math.max(1, values.length - 1) * width, height - 1 - v / max * (height - 2)]);
$: line = points.map(([x, y]) => `${x.toFixed(1)},${y.toFixed(1)}`).join(" ");
$: area = `0,${height} ${line} ${width},${height}`;
</script>

{#if values.length > 1}
	<svg {width} {height} viewBox="0 0 {width} {height}" preserveAspectRatio="none" aria-hidden="true">
		<polygon points={area} class="area" />
		<polyline points={line} class="line" vector-effect="non-scaling-stroke" />
	</svg>
{/if}

<style>
svg {
	display: block;
	overflow: visible;
}
.area {
	fill: var(--bs-primary);
	opacity: 0.15;
}
.line {
	fill: none;
	stroke: var(--bs-primary);
	stroke-width: 1.5;
	stroke-linejoin: round;
}
</style>
