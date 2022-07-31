<script context="module">
import CONFIG from "@/biowasm.json";

export async function load({ params }) {
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	const version = (tool?.versions || []).find(v => v.version === params.version);
	if(!tool)
		return { status: 303, redirect: CONFIG.url };
	if(!version)
		return { status: 303, redirect: `${CONFIG.url}/${tool.name}` };

	return { props: {
		tool,
		version
	} };
}
</script>

<script>
export let tool;
export let version;
</script>

<base href="{CONFIG.url}/{tool.name}/{version.version}/" />

<h3><a href={CONFIG.url}>CDN</a> / <a href="../">{tool.name}</a> / {version.version}</h3>

{#if version.dependencies}
	Dependencies:<br />
	{#each version.dependencies as dependency}
		* <a href="{CONFIG.url}/{dependency.name}/{dependency.version}">{dependency.name} v{dependency.version}</a><br />
	{/each}
	<br />
{/if}

Files:<br />
{#each tool.programs || [tool.name] as program}
	{#each tool.files || ["js", "wasm"] as file}
		<!-- Need full page reload to download the file -->
		* <a sveltekit:reload href="{program}.{file}">{program}.{file}</a><br />
	{/each}
{/each}
