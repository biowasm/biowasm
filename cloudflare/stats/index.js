// This file defines a Cloudflare Worker that:
//   (1) When called as a cron job, tallies up how often each tool was downloaded
//       from the CDN over time, and stores those counts as "summary".
//   (2) When called as a GET endpoint, uses the summary stats calculated by the
//       cron above to plot a chart of download counts over time.

// KV key name
const KEY_SUMMARY = "summary";

// ================================================================
// Web endpoint
// ================================================================
async function handleRequest(request)
{
	// KV format: { <date>: { <tool>: <count> } }
	let stats = JSON.parse(await LOGS.get(KEY_SUMMARY));

	// Convert to format { <tool> : { <date>: <count> } }
	let statsT = {};
	for(let date in stats) {
		for(let tool in stats[date]) {
			const count = stats[date][tool];
			if(!(tool in statsT))
				statsT[tool] = {};
			statsT[tool][date] = count;
		}
	}

	// Convert to Plotly format
	// Plotly format: [ { name: "tool", type: "scatter", x: [date1, date2, ...], y: [count1, count2, ...] } ]
	let data = [];
	for(let prgm in statsT) {
		let x = Object.keys(statsT[prgm]).sort();
		let y = x.map(d => statsT[prgm][d]);
		data.push({
			name: prgm,
			x: x, y: y,
			type: "scatter"
		});
	}

	return new Response(`
        <!doctype html>
        <html lang="en">
          <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
            <meta name="description" content="A repository of C/C++ genomics tools, pre-compiled to WebAssembly for use in a browser">
            <title>biowasm.com</title>
            <link rel="icon" href="favicon.ico">
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.5.3/dist/css/bootstrap.min.css" integrity="sha384-TX8t27EcRE3e/ihU7zmQxVncDAy5uIKz4rEkgIXeMed4M0jlfIDPvg6uqKI2xXr2" crossorigin="anonymous">
          </head>
        
          <body>
            <nav class="navbar navbar-expand-md navbar-dark fixed-top bg-dark">
              <div class="container">
              <a class="navbar-brand" href="https://biowasm.com">biowasm</a>
              <button class="navbar-toggler collapsed" type="button" data-toggle="collapse" data-target="#navbar" aria-controls="navbar" aria-expanded="false" aria-label="Toggle navigation">
                <span class="navbar-toggler-icon"></span>
              </button>
        
              <div id="navbar" class="collapse navbar-collapse">
                <ul class="navbar-nav mr-auto"></ul>
                <ul class="navbar-nav">
                  <li class="nav-item mr-2">
                    <a class="btn btn-primary" href="https://play.biowasm.com" target="_blank" onclick="javascript:window.location=window.location.href.includes('stg') ? 'https://play-stg.biowasm.com' : 'https://play.biowasm.com'; return false">Launch Playground</a>
                  </li>
                  <li class="nav-item mr-2">
                    <a class="nav-link text-light-50" href="https://cdn.biowasm.com" target="_blank" onclick="javascript:window.location=window.location.href.includes('stg') ? 'https://cdn-stg.biowasm.com' : 'https://cdn.biowasm.com'; return false">Packages</a>
                  </li>
                  <li class="nav-item mr-2">
                    <a class="nav-link text-light" href="https://stats.biowasm.com" target="_blank" onclick="javascript:window.location=window.location.href.includes('stg') ? 'https://stats-stg.biowasm.com' : 'https://stats.biowasm.com'; return false">Stats</a>
                  </li>
                  <li class="nav-item">
                    <a class="nav-link text-light-50" href="https://github.com/biowasm" target="_blank">GitHub</a>
                  </li>
                </ul>
              </div>
              </div>
            </nav>
        
            <main role="main">
              <div class="jumbotron mt-5 mt-md-2 pb-1">
                <div class="container mt-4 mb-4">
                  <h2>Stats</h2>
                  <p class="lead">
                    Package downloads over time
                  </p>
                </div>
              </div>
        
              <div class="container">
                <div class="row mt-4">
                  <div class="col-md-12">
  
                    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                    <div id="plot"></div>

                    <script type="application/javascript">
                    var data = ${JSON.stringify(data)};
                    Plotly.newPlot("plot", data);
                    </script>
                  </div>
                </div>
              </div>
            </main>
          </body>
        </html>
    `, { headers: { "content-type": "text/html;charset=UTF-8" } });
}

// ================================================================
// Cron job
// ================================================================
async function handleSchedule(scheduledDate)
{
	console.log(scheduledDate);
	// ========================================================================
	// First, look for raw, unaggregated log events, and aggregate them
	// Raw events: raw:<date>:<tool>/<version>/<program>.wasm:<uuid> value="" metadata={}
	// ========================================================================
	let logsRaw = await getLogs({ prefix: "raw:" });
	let aggregated = aggregateLogs(logsRaw);

	// Now that they're aggregated, we can update the key/value store counts
	for(let date in aggregated) {
		for(let tool in aggregated[date]) {
			const key = `aggregate:${date}:${tool}`;
			const count = parseInt((await LOGS.getWithMetadata(key)).metadata || 0);
			const countToAdd = aggregated[date][tool];
			await LOGS.put(key, "", { metadata: count + countToAdd });
		}
	}

	// And we can delete the raw log events
	for(let log of logsRaw) {
		console.log(`Deleting ${log.name}...`)
		await LOGS.delete(log.name);
	}

	// ========================================================================
	// Now fetch *all* aggregate KVs and save that to a summary KV
	// Aggregate events: aggregate:<date>:<prgm>  value="" metadata=count
	// Note that we store the summary JSON as { "tool": { "date" } } so it's 
	// easier to plot as different series, but we store the key/value pairs as
	// "date:tool" for sorting purposes.
	// ========================================================================
    let logsAggregate = await getLogs({ prefix: "aggregate:" });
    aggregated = aggregateLogs(logsAggregate);

	// Save it so we can fetch it later quickly from the stats page
	await LOGS.put(KEY_SUMMARY, JSON.stringify(aggregated));
}


// ================================================================
// Utility functions
// ================================================================

async function getLogs(config)
{
	// Loop through all keys and handle pagination
	let processed = [];
	let response = { cursor: null };
	do {
		response = await LOGS.list({
			prefix: config.prefix,
			cursor: response.cursor
		});

		processed = processed.concat(response.keys);
	} while(response.cursor != null);

	return processed;
}

function aggregateLogs(logs)
{
	// End result: aggregated = { <date>: { <tool>: <count> } }
	let aggregated = {};
	logs.map(key => {
		// Parse key name: <raw|aggregate>:<date>:<tool>[:uuid]
		const value = key.name.split(":");
		const date = value[1];
		const tool = value[2];
		const countToAdd = key.metadata || 1;

		// Aggregate all versions
		if(!(date in aggregated))
			aggregated[date] = {};
		if(!(tool in aggregated[date]))
			aggregated[date][tool] = 0;
		aggregated[date][tool] += countToAdd;
	});
	return aggregated;
}


// ================================================================
// Event listeners
// ================================================================

addEventListener("fetch", event => {
	return event.respondWith(handleRequest(event.request));
});

addEventListener("scheduled", event => {
	event.waitUntil(handleSchedule(event.scheduledTime));
});
