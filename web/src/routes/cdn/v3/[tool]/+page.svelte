<script context="module">
import { ListGroup, ListGroupItem } from "sveltestrap";
import { getToolURL } from "$lib/utils";
</script>

<script>
export let data;
</script>

<!-- Description -->
<p class="lead">
	{data.tool.description}
</p>

<h5>Versions</h5>

<ListGroup>
	{#each data.tool.versions.reverse() as version}
		<ListGroupItem tag="a" href={getToolURL(data.tool.name, version.version)} action>{version.version}</ListGroupItem>
	{/each}
</ListGroup>
