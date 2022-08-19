<script>
import { page } from "$app/stores";
import { Styles, Container, Collapse, Navbar, NavbarToggler, NavbarBrand, Nav, NavItem, NavLink } from "sveltestrap";
import CONFIG from "@/biowasm.json";

let isOpen = false;
$: path = $page.url.pathname.split("/")[1];
</script>

<!-- Bootstrap CSS and icons -->
<Styles />

<!-- Navbar -->
<Navbar color="light" expand="md" light container>
	<NavbarBrand href="/">
		<img alt="biowasm logo" height="40" src="/logo.png" />
		biowasm
	</NavbarBrand>
	<NavbarToggler on:click={() => isOpen = !isOpen} />
	<Collapse {isOpen} navbar expand="md" on:update={evt => isOpen = evt.detail.isOpen}>
		<Nav class="ms-auto" navbar>
			<NavItem>
				<NavLink href="/" active={path === ""}>Home</NavLink>
			</NavItem>
			<NavItem>
				<NavLink href="/documentation" active={path === "documentation"}>Documentation</NavLink>
			</NavItem>
			<NavItem>
				<NavLink href={CONFIG.url} active={path === "cdn"}>Packages</NavLink>
			</NavItem>
			<NavItem>
				<NavLink href="/stats" active={path === "stats"}>Stats</NavLink>
			</NavItem>
		</Nav>
	</Collapse>
</Navbar>

<!-- Page Content -->
<Container class="mt-4">
	<slot />
</Container>
