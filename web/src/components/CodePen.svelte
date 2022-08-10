<script>
// CodePen IDE - Custom theme ID# 41857:
// * https://codepen.io/embed-theme-builder
// * https://codepen.io/assets/

// Import code samples dynamically!
const codeSamples = import.meta.glob("@/tools/**/examples/*.html", { as: "raw", eager: true });

export let tool;
export let version;
export let lang = "html";
export let code = codeSamples[`../tools/${tool}/examples/${version}.html`];
</script>

<svelte:head>
	<script async src="https://static.codepen.io/assets/embed/ei.js"></script>
</svelte:head>

<div>
	<div class="codepen" data-prefill data-height="350" data-theme-id="41857" data-default-tab="html,result" data-editable="true">
		<pre data-lang={lang}>{code}</pre>
	</div>
</div>
