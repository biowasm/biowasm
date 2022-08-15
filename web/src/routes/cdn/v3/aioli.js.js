import CONFIG from "@/biowasm.json";

// Point /aioli.js to latest Aioli version
export async function GET({ request }) {
	const origin = new URL(request.url).origin;
	return await fetch(`${origin}${CONFIG.url}/aioli/latest/aioli.js`);
}
