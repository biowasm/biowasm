// The biowasm CDN is hosted using Cloudflare Worker Sites where CDN files
// are stored in a key-value store (Cloudflare Worker KV). For example:
//     key=samtools/1.10/samtools.wasm --> value=<samtools.wasm contents>
// 
// A Cloudflare Worker (a serverless function) is the entry point for
// retrieving those files from the key-value store. This index.js file
// defines the code for that entry point. The code is mostly using the
// Cloudflare Workers Site template, but modified to enable CORS (so that
// biowasm modules can be loaded from non-biowasm.com domains!) and to
// log basic stats about the number of downloads per module.

import { getAssetFromKV, mapRequestToAsset } from "@cloudflare/kv-asset-handler"

/**
 * The DEBUG flag will do two things that help during development:
 * 1. we will skip caching on the edge, which makes it easier to
 *    debug.
 * 2. we will return an error message on exception in your Response rather
 *    than the default 404.html page.
 */
const DEBUG = false;

addEventListener("fetch", event => {
  let response = {};

  // Process user request
  try {
    response = handleEvent(event);
  } catch (e) {
    if (DEBUG)
      response = new Response(e.message || e.toString(), { status: 500 });
    response = new Response("Internal Error", { status: 500 });
  }

  // Log basic stats about number of times a .js file was requested *after* we return a response to the user.
  // Documentation: https://github.com/Logflare/cloudflare-app/blob/e5bb250b13d3fbad35e3f87cbcb7e32b35984ee6/workers/utils.js#L18
  let url = new URL(event.request.url);
  if(url.host.startsWith("cdn") && url.host.endsWith(".biowasm.com") && url.pathname.endsWith(".js"))
  {
    async function logEvent(path) {
      // ISO Date Format: <YYYY-MM-DDTHH:mm:ss.sssZ>
      let key = `${new Date().toISOString().split("T").shift()}|${path}`;
      // Increment count
      let counter = parseInt((await LOGS.getWithMetadata(key)).metadata || 0) + 1;
      // Save count in the metadata so we can retrieve it in bulk when doing .list() for CDN stats
      await LOGS.put(key, "", { metadata: counter });
    }
    event.waitUntil(logEvent(url.pathname));
  }

  // Return result
  event.respondWith(response);
})

async function handleEvent(event) {
  const url = new URL(event.request.url);
  let options = {};

  /**
   * You can add custom logic to how we fetch your assets
   * by configuring the function `mapRequestToAsset`
   */
  options.mapRequestToAsset = handlePrefix(/^\//);

  try {
    if (DEBUG) {
      // customize caching
      options.cacheControl = {
        bypassCache: true,
      };
    }
    let response = await getAssetFromKV(event, options);

    // Enable CORS
    response.headers.set("Access-Control-Allow-Origin", "*");
    return response;
  } catch (e) {
    // if an error is thrown try to serve the asset at 404.html
    if (!DEBUG) {
      try {
        let notFoundResponse = await getAssetFromKV(event, {
          mapRequestToAsset: req => new Request(`${new URL(req.url).origin}/404.html`, req),
        });

        return new Response(notFoundResponse.body, { ...notFoundResponse, status: 404 });
      } catch (e) {}
    }

    return new Response(e.message || e.toString(), { status: 500 });
  }
}

/**
 * Here's one example of how to modify a request to
 * remove a specific prefix, in this case `/docs` from
 * the url. This can be useful if you are deploying to a
 * route on a zone, or if you only want your static content
 * to exist at a specific path.
 */
function handlePrefix(prefix) {
  return request => {
    // https://github.com/cloudflare/kv-asset-handler/blob/3949c2481190485ab2c5779031e4dae322a155c6/src/index.ts#L14
    const parsedUrl = new URL(request.url);
    let pathname = parsedUrl.pathname;

    // If path looks like a directory append index.html
    // e.g. If path is /about/ -> /about/index.html
    if (pathname.endsWith("/")) {
      pathname = pathname.concat("index.html");
    }

    parsedUrl.pathname = pathname;
    return new Request(parsedUrl.toString(), request);
  }
}
