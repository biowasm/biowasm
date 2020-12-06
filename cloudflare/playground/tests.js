const assert = require("assert");
const webdriver = require("selenium-webdriver");
const by = webdriver.By;
const until = webdriver.until;

// Config
const TOOLS = ["bedtools", "samtools"];
const URL_HOST = "http://localhost:5000?tool=bedtools";
const DELAY = 300;


// Click "Examples" button and wait till dropdown menu appears
async function showExamples(driver) {
    await driver.findElement(by.xpath("//button[text() = 'Examples']")).click();
    await driver.wait(until.elementLocated(by.className("dropdown-item")));
}

function trimHTML(str) {
    return str.replace(/<([^>]+?)([^>]*?)>(.*?)<\/\1>/ig, "");
}

// Tests
describe("play.biowasm.com", () => {
    let driver;
    let elTools;
    let elCommand;
    let elOutput;
    let tools = {};

    // Setup
    before(async() => {
        // Launch browser
        driver = new webdriver.Builder().forBrowser("chrome").build();
        await driver.get(URL_HOST);

        // Track elements of interest
        elTools = await driver.findElements(by.css(".btn-tool"));
        elCommand = await driver.findElement(by.css("input"));
        elOutput = await driver.findElement(by.css("pre"));

        // Make sure we test all the tools listed on the UI
        for(let tool of elTools) {
            let toolName = await tool.getText();
            tools[toolName] = tool;
            assert(TOOLS.includes(toolName));
        }

        await driver.sleep(DELAY);
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
            await driver.sleep(DELAY);
            await showExamples(driver);

            // Try all the examples for the current tool
            let examples = await driver.findElements(by.className("dropdown-item"));
            for(example of examples)
            {
                // Run example
                await example.click();
                await driver.sleep(DELAY);

                // Get output
                let command = await elCommand.getAttribute("value");
                let output = trimHTML(await elOutput.getAttribute("innerHTML"));
                let outputArr = output.split("\n");

                // Basic validation
                if(command.includes("--help"))
                    assert(output.includes(tool), "Test --help failed");
                else {
                    if(tool == "bedtools")
                        assert(output.includes("chr1") || output.includes("DRR016846"), "Test for bedtools failed");
                    else if(tool == "samtools")
                        assert(output.includes("DRR016846") || output.includes("GL000199") || output.includes("777"), "Test for samtools failed");
                    else
                        assert(false, `Missing tests for ${tool}`);
                }

                // Test piping by re-running command with "| head"
                await elCommand.clear();
                await elCommand.sendKeys(`${command} | head`, webdriver.Key.ENTER);
                await driver.sleep(DELAY);
                // Make sure we get the number of lines we expected
                let outputHead = trimHTML(await elOutput.getAttribute("innerHTML"));
                let outputHeadArr = outputHead.split("\n");
                assert(outputHead == outputArr.slice(0, 10).join("\n"), "Piping to head gave wrong output");
                assert.equal(outputHeadArr.length, Math.min(outputArr.length, 10), "Piping to head gave wrong number of lines");

                // Test file redirection by re-running command with "| head > test"
                await elCommand.clear();
                await elCommand.sendKeys(`${command} | head > test`, webdriver.Key.ENTER);
                await driver.sleep(DELAY);
                await elCommand.clear();
                await elCommand.sendKeys(`cat /data/test`, webdriver.Key.ENTER);
                await driver.sleep(DELAY);
                // Expect the same output
                let outputHeadFile = trimHTML(await elOutput.getAttribute("innerHTML"));
                assert.equal(outputHeadFile, outputHead, "Redirection gave wrong output");

                await showExamples(driver);
            }

            await example.sendKeys(webdriver.Key.ESCAPE);
        });
    }
});
