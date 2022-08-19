<script context="module">
import CONFIG from "@/biowasm.json";
import { ListGroup, ListGroupItem } from "sveltestrap";
</script>

<script>
export let data;
</script>

<base href="{CONFIG.url}/{data.tool.name}/" />

<!-- Description -->
<p class="lead">
	{data.tool.description}
</p>

<h5>Versions</h5>

<ListGroup>
	{#each data.tool.versions.reverse() as version}
		<ListGroupItem tag="a" href="{version.version}" action>{version.version}</ListGroupItem>
	{/each}
</ListGroup>
