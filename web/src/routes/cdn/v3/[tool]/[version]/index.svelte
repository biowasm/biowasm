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
export let tool;
export let version;
</script>

<base href="{CONFIG.url}/{tool.name}/{version.version}/" />

<!-- Dependencies -->
{#each version.dependencies || [] as dependency}
	<a href="{CONFIG.url}/{dependency.name}/{dependency.version}">
		<Badge pill color="primary" class="mb-4">Depends on {dependency.name} v{dependency.version}</Badge>
	</a>
{/each}

<!-- Sample Usage -->
<h5>Sample Usage</h5>


<!-- Files to download -->
<h5>Files</h5>
{#each tool.programs || [tool.name] as program}
	<h6 class="mt-3">{program}</h6>
	<ListGroup>
		{#each tool.files || ["js", "wasm"] as extension}
			<!-- Need full page reload to download the file -->
			<ListGroupItem tag="a" href="{program}.{extension}" action>{program}.{extension}</ListGroupItem>
		{/each}
	</ListGroup>
{/each}
