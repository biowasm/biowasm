import { sveltekit } from "@sveltejs/kit/vite";
import path from "path";

const config = {
	plugins: [sveltekit()],
	resolve: {
		alias: {
			"@": path.resolve(__dirname, "../")
		}
	}
};

export default config;
