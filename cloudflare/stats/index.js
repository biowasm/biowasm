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
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <div id="plot"></div>

        <script type="application/javascript">
        var data = ${JSON.stringify(data)};
        Plotly.newPlot("plot", data);
        </script>  
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
