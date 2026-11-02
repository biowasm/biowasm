<script>
import { browser } from "$app/environment";
import { page } from "$app/stores";
import { Badge, Breadcrumb, BreadcrumbItem, Button, Icon } from "@sveltestrap/sveltestrap";
import CONFIG from "@/biowasm.json";

// Parse tool and version info from URL path
$: [tool, version] = $page.url.pathname.split("/").slice(3);
$: title = `Packages${!tool ? "" : ": " + tool}${!version ? "" : " / " + version}`;
$: if(browser) getDownloadStats(tool, version);
let nbDownloads = null;

$: hasSimd = (CONFIG.tools.find(t => t.name === tool)?.programs || []).some(p => p.endsWith("-simd"));

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

async function getDownloadStats(tool, version) {
	if(!tool || !version)
		return;

	const data = await (await fetch(`/api/v3/stats/${tool}/${version}`)).json();
	nbDownloads = data?.stats?.[tool]?.[version].total || 0;
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
				{version}

				<!-- Show stats -->
				{#if tool && version && nbDownloads !== null}
					<span class="text-small align-middle">
						<Badge pill color="secondary" class="ms-1">
							{nbDownloads.toLocaleString()} downloads
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
