<script>
import { Styles, Container, Collapse, Navbar, NavbarToggler, NavbarBrand, Nav, NavItem, NavLink } from "sveltestrap";

let isOpen = false;
</script>

<!-- Bootstrap CSS and icons -->
<Styles />

<!-- Navbar -->
<Navbar color="light" expand="md" light container>
	<NavbarBrand href="/">
		<img alt="biowasm logo" height="40" src="/logo.png" />
		biowasm
	</NavbarBrand>
	<NavbarToggler on:click={() => (isOpen = !isOpen)} />
	<Collapse {isOpen} navbar expand="md" on:update={evt => isOpen = evt.detail.isOpen}>
		<Nav class="ms-auto" navbar>
			<NavItem>
				<NavLink href="/cdn/v3">Packages</NavLink>
			</NavItem>
			<NavItem>
				<NavLink href="/stats">Stats</NavLink>
			</NavItem>
			<NavItem>
				<NavLink href="https://sandbox.bio/playground">Playground</NavLink>
			</NavItem>
		</Nav>
	</Collapse>
</Navbar>

<!-- Page Content -->
<Container class="mt-3">
	<slot />
</Container>
