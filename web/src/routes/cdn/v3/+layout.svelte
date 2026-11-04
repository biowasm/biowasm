<script>
import { goto } from "$app/navigation";
import { page } from "$app/stores";
import { Badge, Breadcrumb, BreadcrumbItem, Button, Icon } from "@sveltestrap/sveltestrap";
import CONFIG from "@/biowasm.json";

export let data;

$: ({ tool, version } = $page.params);
$: title = `Packages${!tool ? "" : ": " + tool}${!version ? "" : " / " + version}`;

// Newest first, for the version picker
$: versions = [...(CONFIG.tools.find(t => t.name === tool)?.versions || [])].reverse().map(v => v.version);

$: hasSimd = (CONFIG.tools.find(t => t.name === tool)?.programs || []).some(p => p.endsWith("-simd"));

// Count for the index page; same set the table lists
const nbPackages = CONFIG.tools.filter(t => t.listed !== false).length;

$: repo = getRepo(tool, version);

// Where to browse the source this version was built from. Each host spells its "browse a
// tag" URL differently, and Savannah isn't browsable at its git:// clone URL at all.
function getRepo(tool, version) {
	const config = CONFIG.tools.find(t => t.name === tool);
	const branch = (config?.versions || []).find(v => v.version === version)?.branch;
	if(!config?.url || !branch)
		return null;

	const base = config.url.replace(/\.git$/, "");
	const { hostname } = new URL(config.url);

	if(hostname === "github.com")
		return { href: `${base}/tree/${branch}`, icon: "github", host: "GitHub" };
	if(hostname.includes("gitlab"))
		return { href: `${base}/-/tree/${branch}`, icon: "gitlab", host: "GitLab" };
	if(hostname === "git.savannah.gnu.org")
		return { href: `https://git.savannah.gnu.org/cgit/${base.split("/").pop()}.git/tree/?h=${branch}`, icon: "git", host: "Savannah" };
	return { href: base, icon: "git", host: hostname };
}
</script>

<svelte:head>
	<title>{title} | Biowasm</title>
</svelte:head>

<div class="d-flex justify-content-between align-items-start">
<h4>
	<Breadcrumb>
		<BreadcrumbItem>
			{#if tool == null}
				Packages
				<small class="fs-6 fw-normal text-muted ms-2">{nbPackages} packages</small>
			{:else}
				<a href="{CONFIG.url}">Packages</a>
			{/if}
		</BreadcrumbItem>

		{#if tool}
			<BreadcrumbItem>
				{#if version == null}
					{tool}
				{:else}
					<a href="{CONFIG.url}/{tool}">{tool}</a>
				{/if}

				{#if hasSimd}
					<span class="text-small align-middle">
						<Badge pill color="success opacity-75" class="ms-1 fw-normal" title="Ships a SIMD-accelerated build, used automatically where supported">simd</Badge>
					</span>
				{/if}
			</BreadcrumbItem>
		{/if}

		{#if version}
			<BreadcrumbItem>
				{#if versions.length > 1}
					<select class="form-select form-select-sm d-inline-block w-auto align-middle" aria-label="Version"
						value={version} on:change={evt => goto(`${CONFIG.url}/${tool}/${evt.target.value}`)}>
						{#each versions as v, i}
							<option value={v}>{v}{i === 0 ? " (latest)" : ""}</option>
						{/each}
					</select>
				{:else}
					{version}
				{/if}

				<!-- Show stats -->
				{#if data.nbDownloads !== null}
					<span class="text-small align-middle">
						<Badge pill color="secondary" class="ms-1">
							{data.nbDownloads.toLocaleString()} downloads
						</Badge>
					</span>
				{/if}
			</BreadcrumbItem>
		{/if}
	</Breadcrumb>
</h4>

<!-- Link to the exact tag or commit this version was built from -->
{#if repo}
	<Button outline color="primary" size="sm" href={repo.href} target="_blank" rel="noopener noreferrer">
		<Icon name={repo.icon} /> Open in {repo.host}
	</Button>
{/if}
</div>

<slot />

<style>
.text-small {
	font-size: 0.6em;
}
</style>
