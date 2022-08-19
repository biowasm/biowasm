import { redirect } from "@sveltejs/kit";
import CONFIG from "@/biowasm.json";

export async function load() {
	throw redirect(301, CONFIG.url);
}
