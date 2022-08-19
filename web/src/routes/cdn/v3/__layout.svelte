<script>
import { browser } from "$app/env";
import { page } from "$app/stores";
import { Badge, Breadcrumb, BreadcrumbItem } from "sveltestrap";
import CONFIG from "@/biowasm.json";

// Parse tool and version info from URL path
$: [tool, version] = $page.url.pathname.split("/").slice(3);
$: title = `Packages${!tool ? "" : ": " + tool}${!version ? "" : " / " + version}`;
$: if(browser) getDownloadStats(tool, version);
let nbDownloads = null;

async function getDownloadStats(tool, version) {
	if(!tool || !version)
		return;

	const data = await (await fetch(`/api/v3/stats/${tool}/${version}`)).json();
	nbDownloads = data?.stats?.[tool]?.[version].total || 0;
}
</script>

<svelte:head>
	<title>{title}</title>
</svelte:head>

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

<slot />

<style>
.text-small {
	font-size: 0.6em;
}
</style>
