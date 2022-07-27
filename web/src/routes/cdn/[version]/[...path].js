export async function GET({ params }) {
	// Redirect fetch request to existing CDNs?
	return {
		status: 303,
		headers: {
			location: `https://cdn.biowasm.com/${params.version}/${params.path}`
		}
	};
}
