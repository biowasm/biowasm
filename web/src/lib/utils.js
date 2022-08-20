import CONFIG from "@/biowasm.json";

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
