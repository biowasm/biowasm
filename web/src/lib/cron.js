// Cron that gathers all stats from Durable Objects and stores a summary in a
// KV pair. Note that BIOWASM_CONFIG is defined in `bin/postbuild.sh`.
entry_default.scheduled = async (event, env, ctx) => {
	console.log("CONFIG1", BIOWASM_CONFIG);
	console.log("event", event);
	console.log("env", env);
	console.log("ctx", ctx);
};
