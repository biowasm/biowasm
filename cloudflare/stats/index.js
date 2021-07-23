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
		if(prgm == "base")
			continue;
		let x = Object.keys(statsT[prgm]).sort();
		let y = x.map(d => statsT[prgm][d]);
		data.push({
			name: prgm,
			x: x, y: y,
			type: "scatter"
		});
	}

	// Fill gaps in the data with 0 so the plot doesn't look weird
	data = data.map(d => {
		const x = d.x;
		const y = d.y;

		let xFinal = [];
		let yFinal = [];

		for(let i in x)
		{
			// Calculate how many days are missing between the current date and the previous one
			const datePrev = new Date(x[i - 1]);
			const dateCurr = new Date(x[i]);
			const dateDiff = (dateCurr - datePrev) / 86400000;  // 86400000 = 1 day
			if(i == 0) {
				xFinal.push(dateCurr.toISOString().slice(0, 10));
				yFinal.push(y[i]);
				continue;
			}

			// Add those dates in
			for(let j = 0; j < dateDiff; j++)
			{
				datePrev.setDate( datePrev.getDate() + 1 );
				const dateToAdd = datePrev.toISOString().slice(0, 10);
				xFinal.push(dateToAdd);
				yFinal.push(j == dateDiff - 1 ? y[i] : 0);
			}
		}

		// d.mode= "markers+lines";
		d.x = xFinal;
		d.y = yFinal;
		return d;
	})

	return new Response(JSON.stringify(data), { status: 200 });
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
