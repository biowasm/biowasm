// Script to migrate from storing the download counts in the "value" field and into
// the "metadata" field instead. This is because "metadata" is included in the response
// to a .list() command, whereas "value" requires a separate .get() request for each! i.e.
// this is much much faster for calculating stats
// 
// Before:
//   { key: "path", value: 123 }
// After:
//   { key: "path", value: "", metadata: 123 }

addEventListener("fetch", event => {
  event.respondWith(handleRequest(event.request));
});

async function handleRequest(request)
{
    let response = { cursor: null };
    do {
        let response = { cursor: null };
            response = await LOGS_STG.list({
            cursor: response.cursor
        });

        let i = 0;
        for(let obj of response.keys) {
            if(i++ % 50 == 0)
                console.log(`${i} / ${response.keys.length}...`)

            // Need to put one metadata at a time. Currently, "wrangler kv:bulk" doesn't
            // support a PUT operation that includes metadata
            if(obj.metadata != null)
                continue;
            const count = await LOGS_STG.get(obj.name);
            await LOGS_STG.put(obj.name, "", { metadata: parseInt(count) })
        }
    } while(response.cursor != null);

    return new Response('done', {status: 200});
}
