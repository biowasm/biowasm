const webdriver = require("selenium-webdriver");
const driver = new webdriver.Builder().forBrowser("chrome").build();
const by = webdriver.By;
const until = webdriver.until;

// Run all tests
async function run()
{
    await driver.get("http://localhost:5000");

    // Test out each tool (e.g. bedtools, samtools)
    let tools = await driver.findElements(by.css(".jumbotron button"));
    for(tool of tools)
    {
        // Choose tool and click on Examples dropdown
        await driver.sleep(200);
        await tool.click();
        await driver.sleep(200);
        await showExamples();

        // Try all the examples for the current tool
        let examples = await driver.findElements(by.className("dropdown-item"));
        for(example of examples)
        {
            await driver.sleep(200);
            await example.click()
            await driver.sleep(200);
            await showExamples();
        }

        await example.sendKeys(webdriver.Key.ESCAPE);
    }

    await driver.quit();
};

// Click "Examples" button and wait till dropdown menu appears
async function showExamples()
{
    await driver.findElement(by.xpath("//button[text() = 'Examples']")).click();
    await driver.wait(until.elementLocated(by.className("dropdown-item")));
}

run();
