import { GET as getCDN } from "../[tool]/[version]/[file]/+server";

// Point /aioli.js to latest Aioli version
export async function GET({ request, platform }) {
	return getCDN({ request, platform, params: {
		tool: "aioli",
		version: "latest",
		file: "aioli.js"
	} });
}
