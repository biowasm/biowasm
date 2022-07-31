<script context="module">
import CONFIG from "@/biowasm.json";

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

<h3><a href={CONFIG.url}>CDN</a> / {tool.name}</h3>

{#each tool.versions as version}
	* <a href={version.version}>{version.version}</a><br />
{/each}
