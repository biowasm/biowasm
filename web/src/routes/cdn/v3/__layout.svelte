<script>
import CONFIG from "@/biowasm.json";
import { page } from "$app/stores";
import { Breadcrumb, BreadcrumbItem } from "sveltestrap";

// Parse tool and version info from URL path
$: [tool, version] = $page.url.pathname.split("/").slice(3);
$: title = `Packages${!tool ? "" : ": " + tool}${!version ? "" : " / " + version}`;
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
			</BreadcrumbItem>
		{/if}
	</Breadcrumb>
</h4>

<slot />
