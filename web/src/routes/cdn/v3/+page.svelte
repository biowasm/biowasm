<script context="module">
import { goto } from "$app/navigation";
import { Table } from "sveltestrap";
import { getToolURL } from "$lib/utils";
</script>

<script>
export let data = { tools: [] };
</script>

<Table hover>
	<thead>
		<tr>
			<th width="20%">Name</th>
			<th>Description</th>
			<th width="20%">Latest Version</th>
		</tr>
	</thead>
	<tbody>
		{#each data.tools.filter(d => d.listed !== false) as t}
			{@const toolURL = t.versions.length === 1 ? getToolURL(t.name, t.versions[0].version) : getToolURL(t.name)}
			<tr on:click={() => goto(toolURL)}>
				<td class="text-primary fw-bold"><a class="text-decoration-none" href={toolURL}>{t.name}</a></td>
				<td>{t.description}</td>
				<td>{t.versions[0].version}</td>
			</tr>
		{/each}
	</tbody>
</Table>

<style>
	td { cursor: pointer; }
</style>
