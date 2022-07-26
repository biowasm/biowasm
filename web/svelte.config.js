import adapter from "@sveltejs/adapter-cloudflare-workers";

const config = {
	kit: {
		adapter: adapter()
	}
};

export default config;
