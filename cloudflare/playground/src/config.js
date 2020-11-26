// Constants
export const TOOLS = {
	"samtools": {
		aioli: { module: "samtools", version: "1.10" },
		queries: [
            {
                header: "Filtering",
                items: [
                    {
                        label: "samtools view",
                        command: "samtools view /urls/A549_CDKN2A.bam 9:21,967,752-21,995,324",
                        tooltip: "Get reads from chr9",
                        description: "Get reads from a genomic region of interest"
                    },
                    {
                        label: "samtools view -q 10",
                        command: "samtools view -q 10 /urls/A549_CDKN2A.bam 9:21,967,752-21,995,324",
                        tooltip: "Filter by mapping quality",
                        description: "Filter out reads that have a mapping quality < 10 (mapping quality = 5th column)"
                    },
                    {
                        label: "samtools view -F 2",
                        command: "samtools view -F 2 /urls/A549_CDKN2A.bam 9:21,967,752-21,995,324",
                        tooltip: "Filter by SAM flags",
                        description: "Filter out reads that don't satisfy a certain flag condition (flag = 2nd column)"
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
                        command: "samtools view -c /urls/A549_CDKN2A.bam 9:21,967,752-21,995,324",
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

