// CDN
const TESTING = true;
export const BIOWASM_URL = `https://cdn${TESTING ? "-stg" : ""}.biowasm.com`;

// Constants
const BAM_FILE = "A549.bam";
const BAI_FILE = "A549.bam.bai";
const BAM_REGION = "9:22,000,000-22,100,000";
const BAM_HOST_URL = window.location.origin;

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
                        label: "samtools view -?",
                        command: "samtools view -?",
                        tooltip: '"view" documentation',
                        description: "Extensive documentation for <code>samtools view</code> parameters"
                    }
                ]
            }
        ],
		files: [
            {
                name: BAM_FILE,
                url: `${BAM_HOST_URL}/data/A549_CDKN2A.bam`
            },
            {
                name: BAI_FILE,
                url: `${BAM_HOST_URL}/data/A549_CDKN2A.bam.bai`
            }
		]
	},
	"bedtools2": {
		aioli: {
            module: "bedtools2",
            version: "2.29.2",
        },
    },
	"bowtie2": {
		aioli: {
            module: "bowtie2",
            version: "2.4.2",
        },
	}
};

// Simple JavaScript utility functions that simulate Aioli objects
export const UTILITIES = {
    // List files from a folder on the virtual file system
    "ls": {
        "exec": (aioli, args) => fs(aioli, "readdir", args)
    },
    //

}

async function fs(aioli, cmd, args)
{
    return {
        "stdout": (await aioli.fs(cmd, args)).join("\n"),
        "stderr": ""
    };
}
