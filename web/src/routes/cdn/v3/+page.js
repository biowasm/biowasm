import CONFIG from "@/biowasm.json";

// `data` holds the download totals from +page.server.js
export async function load({ data }) {
	return {
		...data,
		tools: CONFIG.tools
	};
}
