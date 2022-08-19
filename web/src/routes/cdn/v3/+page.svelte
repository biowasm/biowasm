<script context="module">
import { goto } from "$app/navigation";
import CONFIG from "@/biowasm.json";
import { Table } from "sveltestrap";
</script>

<script>
export let data = { tools: [] };
</script>

<base href="{CONFIG.url}/" />

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
			<tr on:click={() => goto(t.versions.length === 1 ? `${t.name}/${t.versions[0].version}` : t.name)}>
				<td class="text-primary fw-bold">{t.name}</td>
				<td>{t.description}</td>
				<td>{t.versions[0].version}</td>
			</tr>
		{/each}
	</tbody>
</Table>

<style>
	td { cursor: pointer; }
</style>
