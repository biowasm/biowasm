<script context="module">
import { goto } from "$app/navigation";
import CONFIG from "@/biowasm.json";
import { Table } from "sveltestrap";

export async function load() {
	return { props: {
		tools: CONFIG.tools
	} };
}
</script>

<script>
export let tools;
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
		{#each tools.filter(d => d.listed !== false) as t}
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
