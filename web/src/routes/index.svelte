<script>
import { Icon, Table } from "sveltestrap";
import LinkPackage from "$components/LinkPackage.svelte";
</script>

<svelte:head>
	<title>Biowasm</title>
</svelte:head>

<!-- Jumbotron -->
<div class="p-5 mb-4 bg-light rounded-3">
	<div class="container-fluid py-5">
		<h1 class="display-6 fw-bold col-md-8">Run C/C++ genomics tools in the browser</h1>
		<p class="fs-4 col-md-8">Use bioinformatics command-line tools in your web apps with just a few lines of code. Powered by WebAssembly.</p>
		<a href="/documentation" class="btn btn-primary btn-lg" type="button">
			Get Started <Icon name="arrow-right" />
		</a>
	</div>
</div>

<!-- More info -->
<div class="row align-items-md-stretch">
	<div class="col-md-4">
		<div class="h-100 p-5 text-bg-dark rounded-3">
			<h2 class="fw-bold">Supported Tools</h2>
			<p>
				Biowasm supports many C/C++ bioinformatics tools (<LinkPackage tool="samtools" />, <LinkPackage tool="bedtools" />, <LinkPackage tool="bcftools" />, <LinkPackage tool="bowtie2" />, <LinkPackage tool="minimap2" />, <LinkPackage tool="fastp" />, <LinkPackage tool="wgsim" />, etc.)
				and general data wrangling utilities (<LinkPackage tool="jq" />, <LinkPackage tool="gawk" />, <LinkPackage tool="grep" />, <LinkPackage tool="sed" />, etc.)
			</p>
			<a href="/cdn/v3" class="btn btn-outline-primary" type="button">
				Browse all packages <Icon name="arrow-right" />
			</a>
		</div>
	</div>
	<div class="col-md-8">
		<div class="h-100 p-5 rounded-3">
			<h2 class="fw-bold">Who Uses Biowasm?</h2>
			<Table>
				<thead>
					<tr>
						<th>App</th>
						<th>Why Biowasm?</th>
					</tr>
				</thead>
				<tbody>
					<tr>
						<td><a class="fw-bold" href="https://sandbox.bio">sandbox.bio</a></td>
						<td>Run command-line tools in the browser to power interactive tutorials</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://quinlan-lab.github.io/bedqc/">bedqc</a></td>
						<td>Run <code>bedtools</code> in the browser to validate .bed files</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://genomeribbon.com">Ribbon</a></td>
						<td>Run <code>samtools</code> in the browser to parse, estimate coverage and subsample .bam files</td>
					</tr>
				</tbody>
			</Table>
		</div>
	</div>
</div>
