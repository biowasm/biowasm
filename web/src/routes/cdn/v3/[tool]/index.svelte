<script context="module">
import CONFIG from "@/biowasm.json";
import { ListGroup, ListGroupItem } from "sveltestrap";

export async function load({ params }) {
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	if(!tool)
		return { status: 303, redirect: CONFIG.url };

	return { props: { 
		tool
	} };
}
</script>

<script>
export let tool;
</script>

<base href="{CONFIG.url}/{tool.name}/" />

<!-- Description -->
<p class="lead">
	{tool.description}
</p>

<h5>Versions</h5>

<ListGroup>
	{#each tool.versions.reverse() as version}
		<ListGroupItem tag="a" href="{version.version}" action>{version.version}</ListGroupItem>
	{/each}
</ListGroup>
