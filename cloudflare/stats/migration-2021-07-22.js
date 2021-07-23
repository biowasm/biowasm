// Script to copy over the stats from the key-value store for CDN v1 to the key-value store for CDN v2

import { getLogs } from "./index.js";

addEventListener("fetch", event => {
	event.respondWith(handleRequest(event.request));
});

async function handleRequest(request)
{
	let logsRaw = await getLogs({ prefix: "20" }, LOGS_STG_V1);
	for(let log of logsRaw) {
		const name = log.name.split("|");
		const count = log.metadata;

		const date = name[0];
		const tool = name[1].replace("/v2", "").split("/")[1];

		const newKey = `aggregate:${date}:${tool}`;

		console.log(log.name, "\t", `key=${newKey}`, `value=${count}`)

		// Uncomment this before running
		// await LOGS_STG_V2.put(newKey, "", { metadata: count })
	}

	return new Response('done', {status: 200});
}
