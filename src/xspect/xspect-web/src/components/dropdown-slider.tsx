"use client"

import * as React from "react"
import { DropdownMenuCheckboxItemProps } from "@radix-ui/react-dropdown-menu"

import { Button } from "@/components/ui/button"
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import { Slider } from "@/components/ui/slider"
import { ChevronDown } from "lucide-react"


interface DropdownMenuSlider {
  triggerButtonText?: string;
  value: number;
  max?: number;
  step?: number;
  disabled?: boolean;
}

export function DropdownMenuSlider({
  triggerButtonText,
  value,
  onValueChange,
  max,
  step,
  disabled = false,
}: DropdownMenuSlider) {
  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="ghost" disabled={disabled}>{triggerButtonText}<ChevronDown /></Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent className="w-56">
        <Slider value={[value]} onValueChange={onValueChange} max={max} step={step} className="my-3" />
      </DropdownMenuContent>
    </DropdownMenu>
  )
}
