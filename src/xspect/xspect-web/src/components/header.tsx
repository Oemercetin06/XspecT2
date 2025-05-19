import {
    NavigationMenu,
    NavigationMenuItem,
    NavigationMenuList,
} from "@/components/ui/navigation-menu"
import { Link, NavLink } from "react-router-dom"


export default function Header() {

    const navItems = [
        { name: "Classify", href: "classify", external: false },
        { name: "Filter", href: "filter", external: false },
        { name: "Models", href: "models", external: false },
        { name: "Documentation", href: "https://bionf.github.io/XspecT2/index.html", external: true },
    ]

    return (
        <header className="bg-white shadow-sm">
            <div className="container mx-auto px-4 py-4">
                <div className="flex items-center justify-between">
                    <div className="flex items-center text-xl font-bold">
                        <Link to="/">
                            XspecT
                        </Link>
                    </div>
                    <NavigationMenu>
                        <NavigationMenuList>
                            {navItems.map((item) => (
                                <NavigationMenuItem key={item.name} className="text-gray-700 hover:text-gray-900 font-medium">
                                    <NavLink
                                        className="block px-4 py-2 rounded-md transition-colors duration-200"
                                        to={item.href}
                                        {...item.external ? { onClick: (e) => { e.preventDefault(); window.open(item.href, "_blank"); } } : {}}
                                    >
                                        {item.name}
                                    </NavLink>
                                </NavigationMenuItem>
                            ))}
                        </NavigationMenuList>
                    </NavigationMenu>
                </div>
            </div>
        </header>
    )
}
