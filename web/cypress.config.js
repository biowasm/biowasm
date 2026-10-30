import { defineConfig } from "cypress";
import fs from "fs";
import path from "path";

// env=local serves tools out of build/, so the tests need to know which Aioli version that is
const AIOLI_LATEST = JSON.parse(fs.readFileSync("../biowasm.json")).tools.find(t => t.name === "aioli").versions.at(-1).version;

export default defineConfig({
	fixturesFolder: false,
	video: false,
	screenshotOnRunFailure: false,
	defaultCommandTimeout: 60000,
	fileServerFolder: "../",
	expose: { aioliLatest: AIOLI_LATEST },
	e2e: {
		specPattern: 'tests//**/*.cy.js',
		supportFile: false,
		// cy.exec() was removed in Cypress 16, so the tests ask Node to do this instead
		setupNodeEvents(on) {
			on("task", {
				// Every example we expect to have a test for, relative to the repo root.
				// The .stg/.local copies a test wrote are not examples of their own.
				listExamples() {
					const examples = [];
					const tools = path.resolve("../tools");
					for(const tool of fs.readdirSync(tools)) {
						const dir = path.join(tools, tool, "examples");
						if(!fs.existsSync(dir))
							continue;
						for(const f of fs.readdirSync(dir))
							if(f.endsWith(".html") && !/\.(stg|local)\.html$/.test(f))
								examples.push(`tools/${tool}/examples/${f}`);
					}
					for(const f of fs.readdirSync("src/examples"))
						if(f.endsWith(".html") && !/\.(stg|local)\.html$/.test(f))
							examples.push(`web/src/examples/${f}`);
					return examples;
				},

				// Remove the temporary copy a test wrote
				deleteFile(file) {
					fs.rmSync(path.resolve("..", file), { force: true });
					return null;
				}
			});
		}
	},
});
