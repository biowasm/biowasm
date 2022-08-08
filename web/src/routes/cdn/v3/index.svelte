<script context="module">
import CONFIG from "@/biowasm.json";

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

{#each tools as t}
	* <a href={t.name}>{t.name}</a> v{t.versions[0].version}<br />
{/each}
