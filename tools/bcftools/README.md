## bcftools.wasm

### Usage

```html
<script src="https://cdn.biowasm.com/v2/aioli/latest/aioli.js"></script>
<script type="module">
let CLI = await new Aioli("bcftools/1.10");

// Filter VCF by QUAL column value, and only output chromosome, position, reference, and alternate allele
let output = await CLI.exec(`bcftools query -i %QUAL>=60 -f %CHROM:%POS\t%REF\t%ALT\n /bcftools/annotate.vcf`)
document.write(`<pre>${output}</pre>`);
</script>
```

### Patch
- Need to reset `opt` variables so that it works properly when call `main()` multiple times
