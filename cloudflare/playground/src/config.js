// Constants
export const TOOLS = {
    // 9:22,000,000-22,100,000, 9:21,967,752-21,995,324
    "samtools": {
		aioli: { module: "samtools", version: "1.10" },
		queries: [
            {
                header: "Filtering",
                items: [
                    {
                        label: "samtools view",
                        command: "samtools view /urls/A549_CDKN2A.bam 9:22,000,000-22,100,000",
                        tooltip: "Get reads from chr9",
                        description: "Get reads from a genomic region of interest"
                    },
                    {
                        label: "samtools view -q 20",
                        command: "samtools view -q 20 /urls/A549_CDKN2A.bam 9:22,000,000-22,100,000",
                        tooltip: "Filter by mapping quality",
                        description: "Filter out reads that have a mapping quality < 20 (mapping quality = 5th column of SAM file)"
                    },
                    {
                        label: "samtools view -F 4",
                        command: "samtools view -F 4 /urls/A549_CDKN2A.bam 9:22,000,000-22,100,000",
                        tooltip: "Filter out unmapped reads",
                        description: "Use the <code>-F 4</code> flag to exclude reads where the SAM flag <code>4</code> is set (flag = 2nd column of SAM file) &mdash; <a href='https://broadinstitute.github.io/picard/explain-flags.html' target='_blank'>SAM flags utility</a>"
                    },
                    {
                        label: "samtools view -f 4",
                        command: "samtools view -f 4 /urls/A549_CDKN2A.bam 9:22,000,000-22,100,000",
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
                        command: "samtools idxstats /urls/A549_CDKN2A.bam",
                        tooltip: "Stats about the .bai file",
                        description: "Retrieve statistics about the <code>.bai</code> index"
                    },
                    {
                        label: "samtools view -H",
                        command: "samtools view -H /urls/A549_CDKN2A.bam",
                        tooltip: "Output BAM file header",
                        description: "Retrieve the header from the <code>.bam</code> file"
                    },
                    {
                        label: "samtools view -c",
                        command: "samtools view -c /urls/A549_CDKN2A.bam 9:22,000,000-22,100,000",
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

		]
	},
	"bedtools2": {
		aioli: { module: "bedtools2", version: "2.29.2" }
	}
};

