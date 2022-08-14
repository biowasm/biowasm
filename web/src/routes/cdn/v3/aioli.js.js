import CONFIG from "@/biowasm.json";

// Redirect /aioli.js to latest Aioli
export async function GET() {
	return {
		status: 301,
		headers: {
			location: `${CONFIG.url}/aioli/latest/aioli.js`
		}
	};
}
