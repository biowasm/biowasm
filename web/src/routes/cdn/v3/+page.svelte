<script context="module">
import { onMount } from "svelte";
import { goto } from "$app/navigation";
import { Table } from "sveltestrap";
import { getToolURL } from "$lib/utils";
</script>

<script>
export let data = { tools: [] };

// State
let sortColumn = "downloads";
let sortDirection = "desc";
let stats = null;

// Sort tools
$: tools = data.tools.filter(d => d.listed !== false).sort((a, b) => {
	let diff = 0;

	if(sortColumn === "downloads") {
		const downloadsA = getNbDownloads(a);
		const downloadsB = getNbDownloads(b);
		diff = downloadsA - downloadsB;
	} else {
		diff = a[sortColumn].localeCompare(b[sortColumn]);
	}

	return diff * (sortDirection === "asc" ? 1 : -1);
});

// Load stats on load
onMount(async () => {
	stats = (await (await fetch("/api/v3/stats")).json()).stats;
});

// Utilities
function getNbDownloads(tool) {
	let sum = 0;
	const statsPerVersion = stats?.[tool.name];
	for(let version in statsPerVersion)
		sum += statsPerVersion[version]?.total || 0;
	return sum;
}

function sortBy(col) {
	sortColumn = col;
	if(sortDirection === "asc")
		sortDirection = "desc";
	else
		sortDirection = "asc";
}
</script>

<Table hover>
	<thead>
		<tr>
			<th width="20%" on:click={() => sortBy("name")}>Name</th>
			<th>Description</th>
			<th width="10%">Latest</th>
			<th width="10%" on:click={() => sortBy("downloads")}>Downloads</th>
		</tr>
	</thead>
	<tbody>
		{#each tools as t}
			{@const toolURL = t.versions.length === 1 ? getToolURL(t.name, t.versions.at(-1).version) : getToolURL(t.name)}
			<tr on:click={evt => {
				// Support Cmd+Click on table row
				if(evt.metaKey)
					window.open(toolURL, '_blank');
				else
					goto(toolURL);
			}}>
				<td class="text-primary fw-bold">{t.name}</td>
				<td>{t.description}</td>
				<td>{t.versions.at(-1).version}</td>
				{#key stats}
					<td>{getNbDownloads(t).toLocaleString()}</td>
				{/key}
			</tr>
		{/each}
	</tbody>
</Table>

<style>
	td { cursor: pointer; }
</style>
