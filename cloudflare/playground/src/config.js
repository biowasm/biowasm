// CDN
const TESTING = true;
export const BIOWASM_URL = `https://cdn${TESTING ? "-stg" : ""}.biowasm.com`;

// Constants
const URL_HOST = window.location.origin;
const BAM_FILE = "A549.bam";
const BAI_FILE = "A549.bam.bai";
const BAM_REGION = "9:22,000,000-22,100,000";
const BED_CPG_FILE = "cpg.bed";
const BED_EXONS_FILE = "exons.bed";
const BED_GENOME_FILE = "genome.txt";

// Tools
export const TOOLS = {
    "samtools": {
		aioli: {
            module: "samtools",
            version: "1.10"
        },
		queries: [
            {
                header: "Filtering",
                items: [
                    {
                        label: "samtools view",
                        command: `samtools view /urls/${BAM_FILE} ${BAM_REGION}`,
                        tooltip: "Get reads from chr9",
                        description: "Get reads from a genomic region of interest"
                    },
                    {
                        label: "samtools view -q 20",
                        command: `samtools view -q 20 /urls/${BAM_FILE} ${BAM_REGION}`,
                        tooltip: "Filter by mapping quality",
                        description: "Filter out reads with mapping quality < 20 (mapping quality = 5th column of SAM file)"
                    },
                    {
                        label: "samtools view -F 4",
                        command: `samtools view -F 4 /urls/${BAM_FILE} ${BAM_REGION}`,
                        tooltip: "Filter out unmapped reads",
                        description: "Use the <code>-F 4</code> flag to exclude reads where the SAM flag <code>4</code> is set (flag = 2nd column of SAM file) &mdash; <a href='https://broadinstitute.github.io/picard/explain-flags.html' target='_blank'>SAM flags utility</a>"
                    },
                    {
                        label: "samtools view -f 4",
                        command: `samtools view -f 4 /urls/${BAM_FILE} ${BAM_REGION}`,
                        tooltip: "List unmapped reads",
                        description: "Use the <code>-f 4</code> flag to only include reads where the SAM flag <code>4</code> is set (flag = 2nd column of SAM file) &mdash; <a href='https://broadinstitute.github.io/picard/explain-flags.html' target='_blank'>SAM flags utility</a>"
                    }
                ]
            },
            {
                header: "Stats",
                items: [
                    {
                        label: "samtools idxstats",
                        command: `samtools idxstats /urls/${BAM_FILE}`,
                        tooltip: "Stats about the .bai file",
                        description: "Retrieve statistics about the <code>.bai</code> index"
                    },
                    {
                        label: "samtools view -H",
                        command: `samtools view -H /urls/${BAM_FILE}`,
                        tooltip: "Output BAM file header",
                        description: "Retrieve the header from the <code>.bam</code> file"
                    },
                    {
                        label: "samtools view -c",
                        command: `samtools view -c /urls/${BAM_FILE} ${BAM_REGION}`,
                        tooltip: "Output number of reads",
                        description: "Only returns the number (or <u>c</u>ount) of reads in the region of interest"
                    }
                ]
            },
            {
                header: "Documentation",
                items: [
                    {
                        label: "samtools --help",
                        command: "samtools --help",
                        tooltip: 'Available commands',
                        description: "All available <code>samtools</code> commands"
                    }
                ]
            }
        ],
		files: [
            {
                name: BAM_FILE,
                url: `${URL_HOST}/data/${BAM_FILE}`
            },
            {
                name: BAI_FILE,
                url: `${URL_HOST}/data/${BAI_FILE}`
            }
		]
	},
	"bedtools": {
		aioli: {
            module: "bedtools2",
            version: "2.29.2",
        },
        queries: [
            {
                header: "Intersect",
                items: [
                    {
                        label: "bedtools intersect",
                        command: "bedtools intersect -a cpg.bed -b exons.bed",
                        tooltip: "Get overlapping regions",
                        description: "Find overlapping regions between files <code>cpg.bed</code> and <code>exons.bed</code>"
                    },
                    {
                        label: "bedtools intersect -v",
                        command: "bedtools intersect -a cpg.bed -b exons.bed -v",
                        tooltip: "Get non overlapping regions",
                        description: "The <code>-v</code> flag shows regions in <code>cpg.bed</code> that don't overlap with <code>exons.bed</code>"
                    },
                    {
                        label: "bedtools intersect -wo",
                        command: "bedtools intersect -a cpg.bed -b exons.bed -wo",
                        tooltip: "Calculate amount of overlap",
                        description: "Using the <code>-wo</code> flag, bedtools outputs the amount of overlap (in basepairs) in the last column"
                    },
                    {
                        label: "bedtools intersect -wo -f",
                        command: "bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.7",
                        tooltip: "Require minimum % overlap",
                        description: "The <code>-f</code> flag enforces a minimum overlap for a region to be listed"
                    },
                ]
            },
            {
                header: "Merge",
                items: [
                    {
                        label: "bedtools merge",
                        command: "bedtools merge -i exons.bed",
                        tooltip: "Merge overlapping regions",
                        description: "Output overlapping regions in a .bed file"
                    },
                    {
                        label: "bedtools merge -d",
                        command: "bedtools merge -i exons.bed -d 10e3",
                        tooltip: "Merge nearby regions",
                        description: "Merges regions if they are within 10kb of each other (i.e. not necessarily overlapping)"
                    },
                ]
            },
            {
                header: "Miscellaneous",
                items: [
                    {
                        label: "bedtools complement",
                        command: "bedtools complement -i exons.bed -g genome.txt",
                        tooltip: "Get regions that don't overlap the genome",
                        description: "Get regions of <code>genome.txt</code> that don't overlap with <code>exons.bed</code>. Use <button class='btn btn-sm btn-info terminal' value='cat genome.txt'>cat genome.txt</button> to see the contents of <code>genome.txt</code>."
                    },
                    {
                        label: "bedtools links",
                        command: "bedtools links -i cpg.bed -db hg38",
                        tooltip: "Create links to UCSC",
                        description: "Creates links to the UCSC genome browser for each region"
                    },
                ]
            },
            {
                header: "Documentation",
                items: [
                    {
                        label: "bedtools --help",
                        command: "bedtools --help",
                        tooltip: "Available commands",
                        description: "All available <code>bedtools</code> commands"
                    },
                ]
            },
        ],
		files: [
            {
                name: BED_CPG_FILE,
                url: `${URL_HOST}/data/${BED_CPG_FILE}`
            },
            {
                name: BED_EXONS_FILE,
                url: `${URL_HOST}/data/${BED_EXONS_FILE}`
            },
            {
                name: BED_GENOME_FILE,
                url: `${URL_HOST}/data/${BED_GENOME_FILE}`
            }
		]
    },
	// "bowtie": {
	// 	aioli: {
    //         module: "bowtie2",
    //         version: "2.4.2",
    //     },
    //     queries: [],
    //     files: []
	// }
};

// Simple JavaScript utility functions that can be called from CommandLine.svelte
export const UTILITIES = {
    // Virtual file system utilities
    "pwd"  : async (aioli, args) => fs(aioli, "cwd"),
    "cd"   : async (aioli, args) => fs(aioli, "chdir", args),
    "ls"   : async (aioli, args) => fs(aioli, "readdir", args || await fs(aioli, "cwd")),
    "mv"   : async (aioli, args) => fs(aioli, "rename", ...args.split(" ")),
    "rm"   : async (aioli, args) => fs(aioli, "unlink", args),
    "stat" : async (aioli, args) => fs(aioli, "stat", args),
    "touch": async (aioli, args) => fs(aioli, "writeFile", args, ""),
    "cat"  : async (aioli, args) => fs(aioli, "readFile", args, { encoding: "utf8" }),
    "head" : async (aioli, args) => fs(aioli, "readFile", args, { encoding: "utf8" }).then(d => d.split("\n").slice(0, 10).join("\n")),
    "tail" : async (aioli, args) => fs(aioli, "readFile", args, { encoding: "utf8" }).then(d => d.split("\n").slice(-10).join("\n")),
    // Download a file
    "download": async (aioli, args) => aioli.download(args).then(url => `<strong>${args}:</strong><br />&bullet; <a href="${url}" download=${args}>Download</a><br />&bullet; <a href="${url}" target="_blank">Open in new tab</a> `),
    // Mount a URL
    "mount": async (aioli, args) => aioli.constructor.mount(args, name=args.split("/").pop()).then(() => "ok"),
    // Other utilities
    "echo" : async (aioli, args) => args,
};

async function fs(aioli, cmd, ...args)
{
    let output = await aioli.fs(cmd, ...args);
    if(Array.isArray(output))
        output = output.join("\n");
    else if(typeof output === "object")
        output = JSON.stringify(output, null, 2);
    return output;
}

// Support for simple piping commands
export const PIPING = {
    "tail": d => d.split("\n").slice(0, 10).join("\n"),
    "head": d => d.split("\n").slice(-10).join("\n"),
    "wc -l": d => String(d.trim().split("\n").length)
};
