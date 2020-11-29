const assert = require("assert");
const webdriver = require("selenium-webdriver");
const by = webdriver.By;
const until = webdriver.until;

// Config
const TOOLS = ["bedtools", "samtools"];
const URL_HOST = "http://localhost:5000";
const DELAY = 300;


// Click "Examples" button and wait till dropdown menu appears
async function showExamples(driver)
{
    await driver.findElement(by.xpath("//button[text() = 'Examples']")).click();
    await driver.wait(until.elementLocated(by.className("dropdown-item")));
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
        elTools = await driver.findElements(by.css(".jumbotron button"));
        elCommand = await driver.findElement(by.css("input"));
        elOutput = await driver.findElement(by.css("pre"));

        // Make sure we test all the tools listed on the UI
        for(let tool of elTools) {
            let toolName = await tool.getText();
            tools[toolName] = tool;
            assert(TOOLS.includes(toolName));
        }
        
        await driver.sleep(DELAY * 2);
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
                let output = await elOutput.getAttribute("innerHTML");
                let outputArr = output.split("\n");

                // Basic validation
                if(command.includes("--help"))
                    assert(output.includes(tool));
                else {
                    if(tool == "bedtools")
                        assert(output.includes("chr1"));
                    else if(tool == "samtools")
                        assert(output.includes("DRR016846") || output.includes("GL000199") || output.includes("777"));
                    else
                        assert(false, `Missing tests for ${tool}`);
                }

                await showExamples(driver);
            }

            await example.sendKeys(webdriver.Key.ESCAPE);
        });
    }
});
