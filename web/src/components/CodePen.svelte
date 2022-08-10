<script>
export let code = "Hello World";
export let lang = "html";
</script>

<svelte:head>
	<script async src="https://static.codepen.io/assets/embed/ei.js" on:load={d => {
		console.log("LOADED", d)
		// window.__CPEmbed(".codepen-later");
	}}></script>
</svelte:head>

<div>
	<div class="codepen" data-prefill data-height="270" data-theme-id="39912" data-default-tab="html,result" data-editable="true">
		<pre data-lang={lang}>{code}</pre>
	</div>
</div>
