<script context="module">
import CONFIG from "@/biowasm.json";
import { Badge, ListGroup, ListGroupItem } from "sveltestrap";

export async function load({ params }) {
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	const version = (tool?.versions || []).find(v => v.version === params.version);
	if(!tool)
		return { status: 303, redirect: CONFIG.url };
	if(!version)
		return { status: 303, redirect: `${CONFIG.url}/${tool.name}` };

	return { props: {
		tool,
		version
	} };
}
</script>

<script>
import { onMount } from "svelte";

export let tool;
export let version;

// Load stats on page load
let stats = {};
onMount(async () => {
	const response = await fetch(`/api/v3/stats/${tool.name}/${version.version}`).then(d => d.json());
	const statsPerProgram = response.stats[tool.name];  // cat: { version: { date: ..., total: ... } }
	for(let program in statsPerProgram)
		stats[program] = statsPerProgram[program][version.version]?.total || 0;
})
</script>

<base href="{CONFIG.url}/{tool.name}/{version.version}/" />

<!-- Dependencies -->
{#each version.dependencies || [] as dependency}
	<Badge pill color="primary" class="mb-4">
		Depends on
		<a class="text-light" href="{CONFIG.url}/{dependency.name}/{dependency.version}">
			{dependency.name} v{dependency.version}
		</a>
	</Badge>
{/each}

<!-- Sample Usage -->
<h5>Sample Usage</h5>


<!-- Files to download -->
<h5>Files</h5>
{#each tool.programs || [tool.name] as program}
	<h6 class="mt-3">
		{program}
		{#if stats[program]}
			<Badge pill color="secondary" class="ms-1">{stats[program]} downloads</Badge>
		{/if}
	</h6>
	<ListGroup>
		{#each tool.files || ["js", "wasm"] as extension}
			<!-- Need full page reload to download the file -->
			<ListGroupItem tag="a" href="{program}.{extension}" action>{program}.{extension}</ListGroupItem>
		{/each}
	</ListGroup>
{/each}
