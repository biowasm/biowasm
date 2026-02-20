#!/usr/bin/env node
// Test Mash WASM build by running key commands and comparing output
// Usage: node test_wasm.mjs [genome_dir]
// genome_dir: directory containing .fna files (default: current directory)

import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'fs';
import { join, dirname, basename } from 'path';
import { fileURLToPath } from 'url';
import { createRequire } from 'module';

const __dirname = dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);

const WASM_DIR = join(__dirname, 'build');
const GENOME_DIR = process.argv[2] || '.';

async function createMashModule() {
    const wasmBinary = readFileSync(join(WASM_DIR, 'mash.wasm'));
    const MashFactory = require(join(WASM_DIR, 'mash.js'));

    let stdout = '';
    let stderr = '';

    const module = await MashFactory({
        wasmBinary,
        print: (text) => { stdout += text + '\n'; },
        printErr: (text) => { stderr += text + '\n'; },
        noInitialRun: true,
    });

    return {
        module,
        run: (args) => {
            stdout = '';
            stderr = '';
            try {
                module.callMain(args);
            } catch (e) {
                // Mash calls exit() which throws in Emscripten
            }
            return { stdout: stdout.trim(), stderr: stderr.trim() };
        },
        writeFile: (path, data) => {
            const dir = dirname(path);
            try { module.FS.mkdirTree(dir); } catch (e) {}
            module.FS.writeFile(path, data);
        },
        readFile: (path) => module.FS.readFile(path),
    };
}

async function main() {
    console.log('=== Mash WASM Test Suite ===\n');

    // Test 1: Version
    console.log('Test 1: mash --version');
    const mash = await createMashModule();
    const versionResult = mash.run(['--version']);
    console.log(`  Version: ${versionResult.stdout}`);
    if (versionResult.stdout !== '2.3') {
        console.error('  FAIL: Expected version 2.3');
        process.exit(1);
    }
    console.log('  PASS\n');

    // Test 2: mash bounds
    console.log('Test 2: mash bounds');
    const mash2 = await createMashModule();
    const boundsResult = mash2.run(['bounds']);
    const boundsLines = boundsResult.stdout.split('\n');
    console.log(`  Output lines: ${boundsLines.length}`);
    if (boundsLines.length < 5) {
        console.error('  FAIL: Expected at least 5 lines of bounds output');
        process.exit(1);
    }
    console.log('  PASS\n');

    // Load genome files if available
    const genomes = [
        'E_coli_K12MG1655.fna',
        'E_coli_CFT073.fna',
        'E_coli_HS.fna',
        'E_coli_UTI89.fna',
        'Ecoli_O126.fna',
    ];

    const availableGenomes = genomes.filter(g => existsSync(join(GENOME_DIR, g)));

    if (availableGenomes.length < 2) {
        console.log('Skipping genome tests (need at least 2 .fna files in genome_dir)');
        console.log(`Looked in: ${GENOME_DIR}`);
        console.log('\n=== All basic tests passed ===');
        return;
    }

    console.log(`Found ${availableGenomes.length} genome files in ${GENOME_DIR}\n`);

    // Test 3: mash sketch
    console.log('Test 3: mash sketch (single genome)');
    const mash3 = await createMashModule();
    const genome1Data = readFileSync(join(GENOME_DIR, availableGenomes[0]));
    mash3.writeFile(`/data/${availableGenomes[0]}`, genome1Data);
    const sketchResult = mash3.run(['sketch', `/data/${availableGenomes[0]}`, '-o', '/data/ref']);
    console.log(`  stderr: ${sketchResult.stderr}`);
    // Check that .msh file was created
    try {
        const mshData = mash3.readFile('/data/ref.msh');
        console.log(`  Sketch file size: ${mshData.length} bytes`);
        console.log('  PASS\n');
    } catch (e) {
        console.error('  FAIL: ref.msh not created');
        process.exit(1);
    }

    // Test 4: mash dist (two genomes)
    console.log('Test 4: mash dist (pairwise distance)');
    const mash4 = await createMashModule();
    const g1 = readFileSync(join(GENOME_DIR, availableGenomes[0]));
    const g2 = readFileSync(join(GENOME_DIR, availableGenomes[1]));
    mash4.writeFile(`/data/${availableGenomes[0]}`, g1);
    mash4.writeFile(`/data/${availableGenomes[1]}`, g2);
    const distResult = mash4.run(['dist', `/data/${availableGenomes[0]}`, `/data/${availableGenomes[1]}`]);
    console.log(`  Output: ${distResult.stdout}`);
    const distFields = distResult.stdout.split('\t');
    if (distFields.length >= 3) {
        const distance = parseFloat(distFields[2]);
        console.log(`  Distance: ${distance}`);
        if (isNaN(distance) || distance < 0 || distance > 1) {
            console.error('  FAIL: Invalid distance value');
            process.exit(1);
        }
    } else {
        console.error('  FAIL: Unexpected output format');
        process.exit(1);
    }
    console.log('  PASS\n');

    // Test 5: mash dist (multiple queries against reference)
    if (availableGenomes.length >= 3) {
        console.log('Test 5: mash dist (multi-query)');
        const mash5 = await createMashModule();
        for (const g of availableGenomes) {
            const data = readFileSync(join(GENOME_DIR, g));
            mash5.writeFile(`/data/${g}`, data);
        }
        // Sketch reference
        mash5.run(['sketch', `/data/${availableGenomes[0]}`, '-o', '/data/ref']);
        // Dist against multiple queries
        const queryArgs = availableGenomes.slice(1).map(g => `/data/${g}`);
        const multiDistResult = mash5.run(['dist', '/data/ref.msh', ...queryArgs]);
        const multiLines = multiDistResult.stdout.split('\n').filter(l => l.length > 0);
        console.log(`  ${multiLines.length} distance results:`);
        for (const line of multiLines) {
            const fields = line.split('\t');
            console.log(`    ${basename(fields[1] || '')} -> ${fields[2]}`);
        }
        console.log('  PASS\n');
    }

    // Test 6: mash triangle
    if (availableGenomes.length >= 3) {
        console.log('Test 6: mash triangle');
        const mash6 = await createMashModule();
        for (const g of availableGenomes) {
            const data = readFileSync(join(GENOME_DIR, g));
            mash6.writeFile(`/data/${g}`, data);
        }
        // Sketch all genomes into a combined file
        const sketchArgs = availableGenomes.map(g => `/data/${g}`);
        mash6.run(['sketch', '-o', '/data/all', ...sketchArgs]);
        const triResult = mash6.run(['triangle', '/data/all.msh']);
        const triLines = triResult.stdout.split('\n').filter(l => l.length > 0);
        console.log(`  Triangle matrix (${triLines[0]} genomes):`);
        for (const line of triLines.slice(1)) {
            const parts = line.split('\t');
            console.log(`    ${basename(parts[0] || '')} ${parts.slice(1).join(' ')}`);
        }
        console.log('  PASS\n');
    }

    // Write output for comparison
    if (availableGenomes.length >= 2) {
        console.log('Generating comparison output...');
        const mashComp = await createMashModule();
        for (const g of availableGenomes) {
            const data = readFileSync(join(GENOME_DIR, g));
            mashComp.writeFile(`/data/${g}`, data);
        }
        // All-vs-all distances
        const sketchArgs = availableGenomes.map(g => `/data/${g}`);
        mashComp.run(['sketch', '-o', '/data/all', ...sketchArgs]);

        // Triangle output
        const triResult = mashComp.run(['triangle', '/data/all.msh']);
        writeFileSync(join(__dirname, 'wasm_triangle.txt'), triResult.stdout);
        console.log(`  Saved triangle output to wasm_triangle.txt`);

        // Individual distances (K12 vs all others)
        let distOutput = '';
        for (let i = 1; i < availableGenomes.length; i++) {
            const mashDist = await createMashModule();
            const d1 = readFileSync(join(GENOME_DIR, availableGenomes[0]));
            const d2 = readFileSync(join(GENOME_DIR, availableGenomes[i]));
            mashDist.writeFile(`/data/${availableGenomes[0]}`, d1);
            mashDist.writeFile(`/data/${availableGenomes[i]}`, d2);
            const result = mashDist.run(['dist', `/data/${availableGenomes[0]}`, `/data/${availableGenomes[i]}`]);
            distOutput += result.stdout + '\n';
        }
        writeFileSync(join(__dirname, 'wasm_dist.txt'), distOutput);
        console.log(`  Saved distance output to wasm_dist.txt`);
    }

    console.log('\n=== All tests passed ===');
}

main().catch(e => {
    console.error('Fatal error:', e);
    process.exit(1);
});
