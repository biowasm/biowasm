<script>
// Syntax-highlighted block of code
import hljs from "highlight.js/lib/core";
import javascript from "highlight.js/lib/languages/javascript";
import xml from "highlight.js/lib/languages/xml";
import "highlight.js/styles/github.css";

// Register languages
hljs.registerLanguage("javascript", javascript);
hljs.registerLanguage("xml", xml);

// State
export let code = "<strong>Hello</strong> World.";
export let lang = "xml";

// Trim extra whitespace before highlighting the code
const codeTrimmed = code.split("\n").map(d => d.replace("\t\t", "")).join("\n").trim();
const codeHighlighted = hljs.highlight(codeTrimmed, { language: lang }).value;
</script>

<pre>{@html codeHighlighted}</pre>

<style>
pre {
	margin-left: 20px;
}
</style>
