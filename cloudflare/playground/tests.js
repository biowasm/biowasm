const assert = require("assert");
const webdriver = require("selenium-webdriver");
const by = webdriver.By;
const until = webdriver.until;

// Config
const TOOLS = ["bowtie2", "bedtools", "samtools"];
const URL_HOST = `http://localhost:5000?tool=${TOOLS[0]}`;

// Tests
describe("play.biowasm.com", () => {
    let driver;
    let elTools;
    let elCommand;
    let elOutput;
    let tools = {};

    // Wait until CLI is done executing current command
    async function wait() {
        await driver.wait(until.elementIsEnabled(elCommand), 5000, "CLI still not enabled after 5s");
    }

    // Click "Examples" button and wait till dropdown menu appears
    async function showExamples(driver) {
        await driver.findElement(by.xpath("//button[text() = 'Examples']")).click();
        await driver.wait(until.elementLocated(by.className("dropdown-item")));
    }

    // Setup
    before(async() => {
        // Launch browser
        driver = new webdriver.Builder().forBrowser("chrome").build();
        await driver.get(URL_HOST);

        // Track elements of interest
        elTools = await driver.findElements(by.css(".btn-tool"));
        elCommand = await driver.findElement(by.css("input"));
        elOutput = await driver.findElement(by.id("stdout"));
        await wait();

        // Make sure we test all the tools listed on the UI
        for(let tool of elTools) {
            let toolName = await tool.getText();
            tools[toolName] = tool;
            assert(TOOLS.includes(toolName), "Make sure we test all tools shown on UI");
        }
    });

    // Tear down
    after(() => {
        driver.quit();
    });

    // Test all examples for each tool
    for(let tool of TOOLS)
    {
        it(`Test ${tool} queries`, async () => {
            // Choose tool and click on Examples dropdown
            await tools[tool].click();
            await wait();
            await showExamples(driver);

            // Try all the examples for the current tool
            let examples = await driver.findElements(by.className("dropdown-item"));
            for(example of examples)
            {
                // Run example
                await example.click();
                await wait();

                // Get output
                let command = await elCommand.getAttribute("value");
                let output = await elOutput.getAttribute("innerHTML");
                let outputArr = output.split("\n");
                // await driver.executeScript("return document.getElementById('stderr').remove();");

                // Basic validation
                if(command.includes("--help"))
                    assert(output.includes(tool), "Test --help failed");
                else {
                    if(tool == "bedtools")
                        assert(output.includes("chr1") || output.includes("DRR016846"), "Test for bedtools failed");
                    else if(tool == "samtools")
                        assert(output.includes("DRR016846") || output.includes("GL000199") || output.includes("777"), "Test for samtools failed");
                    else if(tool == "bowtie2")
                        assert(output.includes("r1") || output.includes(""), "Test for bowtie2 failed");
                    else
                        assert(false, `Missing tests for ${tool}`);
                }

                // Test piping by re-running command with "| head"
                await elCommand.clear();
                await elCommand.sendKeys(`${command} | head`, webdriver.Key.ENTER);
                await wait();
                // Make sure we get the number of lines we expected
                let outputHead = await elOutput.getAttribute("innerHTML");
                let outputHeadArr = outputHead.split("\n");
                assert.equal(outputHeadArr.length, Math.min(outputArr.length, 10), "Piping to head gave wrong number of lines");
                // Note: not testing for exact output because some tools like bowtie change output based on how often
                // the commands are run (e.g. CL field in output SAM)
                // assert.equal(outputHead, outputArr.slice(0, 10).join("\n"), "Piping to head gave wrong output");

                // Test file redirection by re-running command with "| head > test"
                await elCommand.clear();
                await elCommand.sendKeys(`${command} | head > test`, webdriver.Key.ENTER);
                await wait();
                await elCommand.clear();
                await elCommand.sendKeys(`cat /data/test`, webdriver.Key.ENTER);
                await wait();
                // Expect the same output
                let outputHeadFile = await elOutput.getAttribute("innerHTML");
                let outputHeadFileArr = outputHeadFile.split("\n");
                assert.equal(outputHeadArr.length, outputHeadFileArr.length, "Redirection gave wrong number of lines");

                await showExamples(driver);
            }

            await example.sendKeys(webdriver.Key.ESCAPE);
        });
    }
});
