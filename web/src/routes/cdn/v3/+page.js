import CONFIG from "@/biowasm.json";

export async function load() {
	return {
		tools: CONFIG.tools
	};
}
