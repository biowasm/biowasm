import adapter from "@sveltejs/adapter-cloudflare";

const config = {
	kit: {
		adapter: adapter()
	}
};

export default config;
