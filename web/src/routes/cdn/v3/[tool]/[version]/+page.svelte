<script context="module">
import { browser } from "$app/env";
import { Badge, Icon, ListGroup, ListGroupItem, Tooltip } from "sveltestrap";
import { getToolURL } from "$lib/utils";
const codeSamples = import.meta.glob("@/tools/**/examples/*.html", { as: "raw", eager: true });  // Import code samples dynamically!
</script>

<script>
import * as ZipJS from "@zip.js/zip.js";
import CodePen from "$components/CodePen.svelte";

export let data;
let busyDownload = false;

$: ({ tool, version, usedBy } = data);
$: code = codeSamples[`../tools/${tool.name}/examples/${version.branch}.html`];

// Download program files as a .zip file
async function downloadAsZip(program) {
	busyDownload = true;

	// Prepare zip file
	const blobWriter = new ZipJS.BlobWriter("application/zip");
	const zipWriter = new ZipJS.ZipWriter(blobWriter);

	// Download and zip up every file
	for(let extension of tool.files || ["js", "wasm"]) {
		const filename = `${program}.${extension}`;
		const url = getToolURL(tool.name, version.version, filename);
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
			<ListGroupItem tag="a" href={getToolURL(t.name, t.version)} action>{t.name} v{t.version}</ListGroupItem>
		{/each}
	</ListGroup>
{/if}

<!-- List dependencies -->
{#if (version.dependencies || []).length > 0}
	<h5 class="mt-4">Depends On</h5>
	<ListGroup>
		{#each version.dependencies as d}
			<ListGroupItem tag="a" href="{getToolURL(d.name, d.version)}" action>{d.name} v{d.version}</ListGroupItem>
		{/each}
	</ListGroup>
{/if}

<!-- Files to download -->
<h5 class="mt-4">
	Files
	<Tooltip target="info-files-{tool.name}">
		You don't need to download the files below if you use the Biowasm CDN. Click for details.
	</Tooltip>
	<!-- Tooltip doesn't work with SSR -->
	{#if browser}
		<a href="/documentation#no-cdn">
			<Icon id="info-files-{tool.name}" name="question-circle-fill" class="text-info" />
		</a>
	{/if}
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
	</h6>
	<ListGroup>
		{#each tool.files || ["js", "wasm"] as extension}
			<ListGroupItem tag="a" href={getToolURL(tool.name, version.version, `${program}.${extension}`)} action>{program}.{extension}</ListGroupItem>
		{/each}
	</ListGroup>
{/each}
