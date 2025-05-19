"use client"
import { PolarAngleAxis, PolarGrid, Radar, RadarChart } from "recharts"
import {
  ChartConfig,
  ChartContainer,
  ChartTooltip,
  ChartTooltipContent,
} from "@/components/ui/chart"

const chartConfig = {
  score: {
    label: "Score",
    color: "hsl(var(--chart-1))",
  },
} satisfies ChartConfig

export interface ResultChartProps {
  data: {
    taxon: string
    score: number
  }[]
}

export function ResultChart({ data }: ResultChartProps) {
  return (
    <ChartContainer
      config={chartConfig}
      className="mx-auto aspect-square max-h-[1000px]"
    >
      <RadarChart data={data}>
        <ChartTooltip cursor={false} content={<ChartTooltipContent />} />
        <PolarAngleAxis dataKey="taxon" />
        <PolarGrid gridType="circle"/>
        <Radar
          dataKey="score"
          fill="var(--primary)"
          fillOpacity={0.6}
          animationDuration={500}
          max={1}
        />
      </RadarChart>
    </ChartContainer>
  )
}
