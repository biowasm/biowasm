import { redirect } from "@sveltejs/kit";
import CONFIG from "@/biowasm.json";

export async function load({ params }) {
	const tool = CONFIG.tools.find(t => t.name === params.tool);
	if(!tool)
		throw redirect(303, CONFIG.url);

	return { 
		tool
	};
}
