// This file defines a Cloudflare Worker that:
//   (1) When called as a cron job, tallies up how often each tool was downloaded
//       from the CDN over time, and stores those counts as "SUMMARY_AGGREGATED".
//   (2) When called as a GET endpoint, uses the summary stats calculated by the
//       cron above to plot a chart of download counts over time.

// KV key names; summaries are calculated and split by versions,
// or aggregateda cross all versions.
const KEY_SUMMARY_SPLIT = "SUMMARY_SPLIT";
const KEY_SUMMARY_AGGREGATED = "SUMMARY_AGGREGATED";

// ================================================================
// Web endpoint
// ================================================================
async function handleRequest(request)
{
  let data = [];
  let stats = JSON.parse(await LOGS.get(KEY_SUMMARY_AGGREGATED));

  // Convert to Plotly format
  // Format: { "seqtk": { "2020-01-01": 4 } }
  for(let prgm in stats)
  {
    let x = Object.keys(stats[prgm]).sort();
    let y = x.map(d => stats[prgm][d]);
    data.push({
      x: x,
      y: y,
      name: prgm,
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
    `, {
    headers: { "content-type": "text/html;charset=UTF-8" }
  });
}

// ================================================================
// Cron job
// ================================================================
async function handleSchedule(scheduledDate)
{
  console.log(scheduledDate);

  let keys = [];             // ["2020-10-01|/aioli/1.3.0/aioli.js", ...]
  let statsSplit = {};       // {"/aioli/1.3.0/aioli.js": { "2020-01-01": 4 }}
  let statsAggregated = {};  // {"aioli": { "2020-01-01": 4 }}

  // Loop through all key names using "cursor" for pagination.
  // Note that the .list() API doesn't return values so we do
  // that separately.
  let response = { cursor: null };
  do {
    // Get next page of keys
    response = await LOGS.list({
      cursor: response.cursor,
      limit: 100
    });

    // Extract key names
    keys = keys.concat(response.keys.map(d => d.name));
  } while(response.cursor != null)

  // Tally up the counts
  for(let key of keys)
  {
    console.log(key);
    let count = + await LOGS.get(key);
    [date, path] = key.split("|");
    if(path == null)
      continue;
    [_, prgm] = path.split("/");

    // Split by version
    if(!(path in statsSplit))
      statsSplit[path] = {};
    statsSplit[path][date] = count;
    // Aggregate all versions
    if(!(prgm in statsAggregated))
      statsAggregated[prgm] = {};
    statsAggregated[prgm][date] = count;
  }

  await LOGS.put(KEY_SUMMARY_SPLIT, JSON.stringify(statsSplit));
  await LOGS.put(KEY_SUMMARY_AGGREGATED, JSON.stringify(statsAggregated));
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
