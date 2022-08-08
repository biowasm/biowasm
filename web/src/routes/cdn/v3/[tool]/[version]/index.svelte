<script context="module">
import CONFIG from "@/biowasm.json";
import { Badge, Icon, ListGroup, ListGroupItem, Tooltip } from "sveltestrap";

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

export let tool;
export let version;
export let usedBy = [];

let stats = {};
let busyDownload = false;

// Load stats on page load
onMount(async () => {
	const response = await fetch(`/api/v3/stats/${tool.name}/${version.version}`).then(d => d.json());
	const statsPerProgram = response.stats[tool.name];  // cat: { version: { date: ..., total: ... } }
	for(let program in statsPerProgram)
		stats[program] = statsPerProgram[program][version.version]?.total || 0;
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

<!-- Sample Usage -->
<h5 class="mt-4">Sample Usage</h5>


{#if usedBy.length > 0}
<h5 class="mt-4">Used By</h5>
<ListGroup>
	{#each usedBy as t}
		<ListGroupItem tag="a" href="{CONFIG.url}/{t.name}/{t.version}" action>{t.name} v{t.version}</ListGroupItem>
	{/each}
</ListGroup>
{/if}

<!-- Files to download -->
<h5 class="mt-4">Files</h5>
{#each tool.programs || [tool.name] as program}
	<h6 class="mt-3">
		{program}
		
		<!-- Button to download all files for this program -->
		<Tooltip target="icon-{program}">
			Download files as .zip file
		</Tooltip>
		<Badge class="ms-2" color={busyDownload ? "secondary" : "primary"} style="cursor: {busyDownload ? "default" : "pointer"}">
			<Icon id="icon-{program}" name="download" onclick={() => downloadAsZip(program)} />
		</Badge>

		<!-- List dependencies -->
		{#each version.dependencies || [] as dependency}
			<Badge pill color="primary" class="ms-1">
				Depends on
				<a class="text-light" href="{CONFIG.url}/{dependency.name}/{dependency.version}">
					{dependency.name} v{dependency.version}
				</a>
			</Badge>
		{/each}

		<!-- Show overall download stats -->
		{#if stats[program]}
			<Badge pill color="secondary" class="ms-1">{stats[program]} downloads</Badge>
		{/if}
	</h6>
	<ListGroup>
		{#each tool.files || ["js", "wasm"] as extension}
			<ListGroupItem tag="a" href="{program}.{extension}" action>{program}.{extension}</ListGroupItem>
		{/each}
	</ListGroup>
{/each}
