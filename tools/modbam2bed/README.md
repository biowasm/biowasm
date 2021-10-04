# modbam2bed.wasm

Summarization of BAM files containing modified base information to BED file with counts.

### Usage

```html
Download sample data from <a href="https://github.com/biowasm/biowasm/files/7261057/biowasm_sample.tar.gz">here</a> and select the files below.<br /><br />
<input id="myfiles" type="file" multiple>

<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
// Note that we use `script type="module"` so we can use top-level await statements
const CLI = await new Aioli("modbam2bed/0.3.1");

async function run(event) {
    // First, mount the file(s) to a virtual file system
    await CLI.mount(event.target.files);

	// Then call modbam2bed
	const output = await CLI.exec("modbam2bed -e -m 5mC chr1_80k.fasta sample.bam -r chr1:60000-70000");
	console.log(output);
}

document.getElementById("myfiles").addEventListener("change", run, false);
</script>
```

### Dependencies
- Depends on a dev version of htslib

### Patch
- Disable threads
