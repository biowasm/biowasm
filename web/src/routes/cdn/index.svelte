<script context="module">
import CONFIG from "@/biowasm.json";

export async function load() {
	return {
		status: 301,
		redirect: CONFIG.url
	};
}
</script>
