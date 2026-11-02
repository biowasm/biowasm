import CONFIG from "@/biowasm.json";
import SIZES from "@/biowasm.sizes.json";

export function getToolURL(toolName, versionName, fileName) {
	let url = `${CONFIG.url}`;
	if(toolName)
		url += `/${toolName}`;
	if(versionName)
		url += `/${versionName}`;
	if(fileName)
		url += `/${fileName}`;
	return url;
}

// Byte size of a published file, recorded by bin/manifest.py at build time
export function getFileSize(toolName, versionName, fileName) {
	return SIZES[`${toolName}/${versionName}/${fileName}`] || null;
}

export function formatSize(bytes) {
	if(!bytes)
		return "";
	return bytes >= 1000 * 1000 ? `${(bytes / 1e6).toFixed(1)} MB` : `${Math.round(bytes / 1000)} KB`;
}
