<script context="module">
import CONFIG from "@/biowasm.json";
import { ListGroup, ListGroupItem } from "sveltestrap";

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

<ListGroup>
	{#each tools as t}
		<ListGroupItem tag="a" href={t.name} action>{t.name} v{t.versions[0].version}</ListGroupItem>
	{/each}
</ListGroup>
