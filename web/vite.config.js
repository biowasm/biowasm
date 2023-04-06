import { sveltekit } from "@sveltejs/kit/vite";
import path from "path";

const config = {
	plugins: [sveltekit()],
	resolve: {
		alias: {
			"@": path.resolve(__dirname, "../"),
			"$components": path.resolve(__dirname, "./src/components"),
			"$examples": path.resolve(__dirname, "./src/examples"),
		}
	},
	ssr: {
		// Avoids "cannot use import statement outside a module" error
		noExternal: ["@popperjs/core"]
	}
};

export default config;
