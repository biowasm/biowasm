<script context="module">
import { onMount } from "svelte";
import { goto } from "$app/navigation";
import { page } from "$app/stores";
import { Badge, Input, Table } from "@sveltestrap/sveltestrap";
import { getToolURL, getFileSize, formatSize } from "$lib/utils";
</script>

<script>
export let data = { tools: [] };

// State
let sortColumn = "downloads";
let sortDirection = "desc";
let stats = null;

// Search and category live in the URL so a filtered view can be linked to
$: query = $page.url.searchParams.get("q") || "";
$: category = $page.url.searchParams.get("category") || "";

$: listed = data.tools.filter(d => d.listed !== false);
$: categories = [...new Set(listed.flatMap(t => t.categories || []))].sort();

// Sort tools
// `stats` is in here on purpose: it arrives after mount, and without it the default
// sort by downloads would run against the zeroes we start with and never rerun
$: tools = listed.filter(t => matchesCategory(t, category)).filter(t => matchesQuery(t, query)).sort((a, b) => {
	let diff = 0;

	if(sortColumn === "downloads") {
		diff = getNbDownloads(a, stats) - getNbDownloads(b, stats);
	} else if(sortColumn === "size") {
		diff = (getSize(a) || 0) - (getSize(b) || 0);
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
// A tool ships SIMD if it publishes a "-simd" variant of one of its programs.
// Aioli picks that variant automatically when the browser supports it.
function hasSimd(tool) {
	return (tool.programs || []).some(p => p.endsWith("-simd"));
}

// Programs a package ships. "-simd" builds aren't listed: Aioli swaps one in for the program
// you asked for, so e.g. bwa ships bwa and bwa-simd but you only ever call bwa
function getPrograms(tool) {
	return (tool.programs || [tool.name]).filter(p => !p.endsWith("-simd"));
}

// What one program of this tool costs to download, at its latest version. Not the sum across
// programs: each is a standalone module, so you only ever fetch the one you asked for
function getSize(tool) {
	const version = tool.versions.at(-1).version;
	const program = getPrograms(tool)[0];
	let sum = 0;
	for(const file of tool.files || ["js", "wasm"])
		sum += getFileSize(tool.name, version, `${program}.${file}`) || 0;
	return sum || null;
}

// Match on program names too, otherwise searching "bgzip" or "hmmsearch" finds nothing even
// though htslib and hmmer ship them
function matchesQuery(tool, query) {
	if(!query)
		return true;
	const q = query.toLowerCase();
	return tool.name.toLowerCase().includes(q)
		|| (tool.description || "").toLowerCase().includes(q)
		|| getPrograms(tool).some(p => p.toLowerCase().includes(q));
}

function matchesCategory(tool, category) {
	return !category || (tool.categories || []).includes(category);
}

function setParam(key, value) {
	const url = new URL($page.url);
	if(value)
		url.searchParams.set(key, value);
	else
		url.searchParams.delete(key);
	goto(url, { replaceState: true, keepFocus: true, noScroll: true });
}

function getNbDownloads(tool, stats) {
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

<div class="d-flex flex-wrap gap-2 align-items-center mb-3">
	<Input
		type="search"
		placeholder="Search tools and descriptions"
		value={query}
		on:input={evt => setParam("q", evt.target.value)}
		class="w-auto flex-grow-1" />

	<div class="d-flex flex-wrap gap-1">
		{#each categories as c}
			<button type="button" class="btn btn-sm {category === c ? 'btn-primary' : 'btn-outline-secondary'}"
				on:click={() => setParam("category", category === c ? "" : c)}>{c}</button>
		{/each}
	</div>
</div>

{#if query || category}
	<p class="text-muted small">
		{tools.length} of {listed.length} packages
		<button type="button" class="btn btn-link btn-sm p-0 ms-1 align-baseline"
			on:click={() => goto($page.url.pathname, { replaceState: true, noScroll: true })}>clear filters</button>
	</p>
{/if}

<Table hover>
	<thead>
		<tr>
			<th width="20%" on:click={() => sortBy("name")}>Name</th>
			<th>Description</th>
			<th width="10%">Latest</th>
			<th width="10%" on:click={() => sortBy("size")}>Size</th>
			<th width="10%" on:click={() => sortBy("downloads")}>Downloads</th>
		</tr>
	</thead>
	<tbody>
		{#each tools as t}
			{@const toolURL = t.versions.length === 1 ? getToolURL(t.name, t.versions.at(-1).version) : getToolURL(t.name)}
			{@const programs = getPrograms(t)}
			<tr on:click={evt => {
				// Support Cmd+Click on table row
				if(evt.metaKey)
					window.open(toolURL, '_blank');
				else
					goto(toolURL);
			}}>
				<td class="text-primary fw-bold">
					{t.name}
					{#if programs.length > 1}
						<span class="text-muted fw-normal" title="{programs.length} programs in this package">({programs.length})</span>
					{/if}
					{#if hasSimd(t)}
						<Badge pill color="success opacity-75" class="ms-1 fw-normal" title="Ships a SIMD-accelerated build, used automatically where supported">simd</Badge>
					{/if}
				</td>
				<td>{t.description}</td>
				<td>{t.versions.at(-1).version}</td>
				<td class="text-muted" title={programs.length > 1 ? "Download size of one program" : ""}>{formatSize(getSize(t))}</td>
				<td>{getNbDownloads(t, stats).toLocaleString()}</td>
			</tr>
		{/each}
	</tbody>
</Table>

<style>
	td { cursor: pointer; }
</style>
