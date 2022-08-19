import { redirect } from "@sveltejs/kit";
import CONFIG from "@/biowasm.json";

export async function load({ params }) {
	// Get tool/version info
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	const version = (tool?.versions || []).find(v => v.version === params.version);
	if(!tool)
		throw redirect(303, CONFIG.url);
	if(!version)
		throw redirect(303, `${CONFIG.url}/${tool.name}`);

	// Find tools that depend on this current tool
	const usedBy = [];
	if(tool && version) {
		// Loop through all tools
		CONFIG.tools.forEach(t => {
			// Loop through each tool's versions
			t.versions.forEach(v => {
				// Are there any that depend on this tool?
				const dependsOnThisTool = (v.dependencies || []).find(d => d.name === tool.name && d.version === version.version);
				if(dependsOnThisTool)
					usedBy.push({
						name: t.name,
						version: v.version
					});
			})
		});
	}

	return {
		tool,
		version,
		usedBy
	};
}
