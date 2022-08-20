import { defineConfig } from "cypress";

export default defineConfig({
	fixturesFolder: false,
	video: false,
	screenshotOnRunFailure: false,
	defaultCommandTimeout: 60000,
	fileServerFolder: "../",
	e2e: {
		experimentalSessionAndOrigin: true,
		specPattern: 'tests//**/*.cy.js',
		supportFile: false
	},
});
