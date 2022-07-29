// Redirect /aioli.js to latest Aioli
export async function GET() {
	return {
		status: 301,
		headers: {
			location: "/cdn/v3/aioli/latest/aioli.js"
		}
	};
}
