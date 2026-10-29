import { defineConfig } from "cypress";
import fs from "fs";
import path from "path";

export default defineConfig({
	fixturesFolder: false,
	video: false,
	screenshotOnRunFailure: false,
	defaultCommandTimeout: 60000,
	fileServerFolder: "../",
	e2e: {
		specPattern: 'tests//**/*.cy.js',
		supportFile: false,
		// cy.exec() was removed in Cypress 16, so the tests ask Node to do this instead
		setupNodeEvents(on) {
			on("task", {
				// Every example we expect to have a test for, relative to the repo root
				listExamples() {
					const examples = [];
					const tools = path.resolve("../tools");
					for(const tool of fs.readdirSync(tools)) {
						const dir = path.join(tools, tool, "examples");
						if(!fs.existsSync(dir))
							continue;
						for(const f of fs.readdirSync(dir))
							if(f.endsWith(".html") && !f.endsWith(".stg.html"))
								examples.push(`tools/${tool}/examples/${f}`);
					}
					for(const f of fs.readdirSync("src/examples"))
						if(f.endsWith(".html") && !f.endsWith(".stg.html"))
							examples.push(`web/src/examples/${f}`);
					return examples;
				},

				// Remove the temporary .stg.html a test wrote
				deleteFile(file) {
					fs.rmSync(path.resolve("..", file), { force: true });
					return null;
				}
			});
		}
	},
});
