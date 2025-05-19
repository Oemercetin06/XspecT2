"use client"
import { DropdownMenuCheckboxItemProps } from "@radix-ui/react-dropdown-menu"

import { Button } from "@/components/ui/button"
import {
  DropdownMenu,
  DropdownMenuCheckboxItem,
  DropdownMenuContent,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import { ChevronDown } from "lucide-react"

type Checked = DropdownMenuCheckboxItemProps["checked"]

interface DropdownCheckboxItem {
  id: string;
  label: string;
  checked: Checked;
  onCheckedChange: (checked: Checked) => void;
  disabled?: boolean;
}

interface DropdownMenuCheckboxesProps {
  triggerButtonText?: string;
  labelText?: string;
  items: DropdownCheckboxItem[];
  disabled?: boolean;
}

export function DropdownMenuCheckboxes({
  triggerButtonText,
  labelText,
  items,
  disabled = false,
}: DropdownMenuCheckboxesProps) {
  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="ghost" disabled={disabled}>{triggerButtonText}<ChevronDown /></Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent className="w-56">
        {labelText && (
          <>
            <DropdownMenuLabel>{labelText}</DropdownMenuLabel>
            <DropdownMenuSeparator />
          </>
        )}
        {items.map((item) => (
          <DropdownMenuCheckboxItem
            key={item.id}
            checked={item.checked}
            onCheckedChange={item.onCheckedChange}
            disabled={item.disabled}
          >
            {item.label}
          </DropdownMenuCheckboxItem>
        ))}
      </DropdownMenuContent>
    </DropdownMenu>
  )
}
