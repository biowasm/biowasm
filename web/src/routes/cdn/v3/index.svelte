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
		{#each tools as t}
			<tr on:click={() => goto(t.name)}>
				<td><a href={t.name}>{t.name}</a></td>
				<td>{t.description}</td>
				<td>v{t.versions[0].version}</td>
			</tr>
		{/each}
	</tbody>
</Table>

<style>
	td { cursor: pointer; }
</style>
