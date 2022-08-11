<script>
import CONFIG from "@/biowasm.json";

export let tool = "seqtk";
export let version = CONFIG.tools.find(t => t.name === tool).versions.at(-1);
</script>

<a class="fw-bold" href="/cdn/v3/{tool}/{version.version}">{tool}</a>
