<script context="module">
import CONFIG from "@/biowasm.json";
import { browser } from "$app/env";
import { Badge, Icon, ListGroup, ListGroupItem, Tooltip } from "sveltestrap";
// Import code samples dynamically!
const codeSamples = import.meta.glob("@/tools/**/examples/*.html", { as: "raw", eager: true });

export async function load({ params }) {
	// Get tool/version info
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	const version = (tool?.versions || []).find(v => v.version === params.version);
	if(!tool)
		return { status: 303, redirect: CONFIG.url };
	if(!version)
		return { status: 303, redirect: `${CONFIG.url}/${tool.name}` };

	// Find tools that depend on this current tool
	const usedBy = [];
	if(tool && version) {
		// Loop through all tools
		CONFIG.tools.forEach(t => {
			// Loop through each tool's versions
			t.versions.forEach(v => {
				// Are there any that depend on this tool?
				const dependsOnThisTool = (v.dependencies || []).find(d => d.name === tool.name && d.version === version.version);
				if(dependsOnThisTool)
					usedBy.push({
						name: t.name,
						version: v.version
					});
			})
		});
	}

	return { props: {
		tool,
		version,
		usedBy
	} };
}
</script>

<script>
import { onMount } from "svelte";
import * as ZipJS from "@zip.js/zip.js";
import { sparkline } from "@fnando/sparkline";
import CodePen from "$components/CodePen.svelte";

export let tool;
export let version;
export let usedBy = [];

let busyDownload = false;
let busyStats = true;
let elSparklines = {};

// Load sample code from this repo
$: code = codeSamples[`../tools/${tool.name}/examples/${version.branch}.html`];

onMount(async () => {
	busyStats = true;

	// Fetch stats for all programs in parallel
	const stats = {};
	const promises = [];
	for(let program of tool.programs) {
		const promise = fetch(`/api/v3/stats/${tool.name}/${version.version}/${program}`)
			.then(d => d.json())
			.then(d => stats[program] = d?.stats?.[tool.name]?.[version.version]?.[program]);
		promises.push(promise);
	}
	await Promise.all(promises);

	// Plot Sparklines
	for(let program of tool.programs) {
		const data = [];
		Object.entries(stats[program]).forEach(xy => {
			if(xy[0] === "total")
				return;
			data.push({
				date: xy[0], value: xy[1]
			});
		});
		console.log(data)
		// const dates = Object.keys(stats[program]).filter(d => d !== "total");
		// console.log(stats[program])
		sparkline(elSparklines[program], data);
	}

	busyStats = false;
});

// Download program files as a .zip file
async function downloadAsZip(program) {
	busyDownload = true;

	// Prepare zip file
	const blobWriter = new ZipJS.BlobWriter("application/zip");
	const zipWriter = new ZipJS.ZipWriter(blobWriter);

	// Download and zip up every file
	for(let extension of tool.files || ["js", "wasm"]) {
		const filename = `${program}.${extension}`;
		const url = `${CONFIG.url}/${tool.name}/${version.version}/${filename}`;
		const blob = await (await fetch(url)).blob();
		await zipWriter.add(filename, new ZipJS.BlobReader(blob));
	}
	await zipWriter.close();

	// Trigger download
	const anchorElement = document.createElement("a");
	anchorElement.href = URL.createObjectURL(blobWriter.blob);
	anchorElement.download = `${program}-${version.version}.zip`;
	anchorElement.click();

	busyDownload = false;
}
</script>

<base href="{CONFIG.url}/{tool.name}/{version.version}/" />

<!-- Description -->
<p class="lead">
	{tool.description}
</p>

<!-- Sample Code (use `browser` check to skip SSR) -->
{#if browser && code}
	<h5 class="mt-4">Sample Code</h5>
	<!-- Re-render CodePen component when change tool -->
	{#key tool}
		<CodePen code={code} />
	{/key}
{/if}

<!-- List who relies on this package -->
{#if usedBy.length > 0}
	<h5 class="mt-4">Used By</h5>
	<ListGroup>
		{#each usedBy as t}
			<ListGroupItem tag="a" href="{CONFIG.url}/{t.name}/{t.version}" action>{t.name} v{t.version}</ListGroupItem>
		{/each}
	</ListGroup>
{/if}

<!-- List dependencies -->
{#if (version.dependencies || []).length > 0}
	<h5 class="mt-4">Depends On</h5>
	<ListGroup>
		{#each version.dependencies as d}
			<ListGroupItem tag="a" href="{CONFIG.url}/{d.name}/{d.version}" action>{d.name} v{d.version}</ListGroupItem>
		{/each}
	</ListGroup>
{/if}

<!-- Files to download -->
<h5 class="mt-4">
	Files
	<Tooltip target="info-files-{tool.name}">
		Only download these files if you don't use the Biowasm CDN. Click for details.
	</Tooltip>
	<a href="/documentation/">
		<Icon id="info-files-{tool.name}" name="question-circle-fill" class="text-info" />
	</a>
</h5>
{#each tool.programs || [tool.name] as program}
	<h6 class="mt-3">
		{program}
		
		<!-- Button to download all files for this program -->
		<span on:click={() => downloadAsZip(program)} >
			<Badge class="ms-2" color={busyDownload ? "secondary" : "primary"} style="cursor: {busyDownload ? "default" : "pointer"}">
				<Icon name="download" />
				<span class="ps-1">Download as .zip</span>
			</Badge>
		</span>

		<!-- Show overall download stats -->
		<!-- {#if !busyStats} -->
			<svg bind:this={elSparklines[program]} class="sparkline" width="100" height="20" stroke-width="3"></svg>
		<!-- {/if} -->
	</h6>
	<ListGroup>
		{#each tool.files || ["js", "wasm"] as extension}
			<ListGroupItem tag="a" href="{program}.{extension}" action>{program}.{extension}</ListGroupItem>
		{/each}
	</ListGroup>
{/each}

<style>
svg {
	stroke: rgba(61, 133, 222);
	fill: rgba(61, 133, 222, 0.3);
}
</style>
