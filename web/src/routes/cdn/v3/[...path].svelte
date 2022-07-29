<script context="module">
import CONFIG from "../../../../../biowasm.json";

// Fetch info before component loads
export async function load({ params }) {
	// Parse URL parameters
	let [tool, version] = params.path.split("/", 2);
	tool ||= null;
	version ||= null;

	return { props: {
		tool, version
	} };
}
</script>

<script>
export let tool;
export let version;

$: tools = CONFIG.tools.filter(t => tool === null || t.name === tool);
// $: dependants = tools.
</script>

<h3>
	<a href="{CONFIG.url}">CDN</a>
	{#if tool}
		/ <a href="{CONFIG.url}/{tool}">{tool}</a>
	{/if}
	{#if version}
		/ {version}
	{/if}
</h3>

{#each tools as t}
	{#if !tool}
		<a href="{CONFIG.url}/{t.name}">
			<strong>{t.name}</strong>
		</a>
		<br />
	{/if}

	<!-- View all versions -->
	{#if version === null}
		{#each t.versions as version}
			* <a href="{CONFIG.url}/{t.name}/{version.version}">{version.version}</a>
			{#if version.dependencies}
				( depends on:
					{#each version.dependencies as dependency}
						<a href="{CONFIG.url}/{dependency.name}/{dependency.version}">{dependency.name} v{dependency.version}</a>
					{/each}
				)
			{/if}
			<br />
		{/each}
		<hr />

	<!-- View one version -->
	{:else}
		{#each t.programs || [tool] as program}
			{#each t.files || ["js", "wasm"] as file}
				<!-- Need full page reload to download the file -->
				* <a sveltekit:reload href="{CONFIG.url}/{t.name}/{version}/{program}.{file}">{program}.{file}</a><br />
			{/each}
		{/each}
	{/if}
{/each}
