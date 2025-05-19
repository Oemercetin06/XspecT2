import { Model, ModelTableEntry } from "@/types/models";
import { ColumnDef } from "@tanstack/react-table"
import { useEffect, useState } from "react";
import { DataTable } from "@/components/data-table"
import { getModels } from "../api";
import { Separator } from "./ui/separator";
import { ArrowRight } from "lucide-react";
import { Link } from "react-router-dom";

export default function Models() {
    const [tableData, setTableData] = useState<ModelTableEntry[] | null>(null);

    useEffect(() => {
        getModels()
            .then((data: Model[]) => {
                const tableEntries: ModelTableEntry[] = [];
                Object.entries(data).forEach(([modelType, modelList]) => {
                    (modelList as string[]).forEach((modelName: string) => {
                        tableEntries.push({ model_type: modelType, name: modelName });
                    });
                });
                setTableData(tableEntries);
            })
            .catch((err) => console.error('Fetch error:', err));
    }, []);

    const columns: ColumnDef<ModelTableEntry>[] = [
        { header: "Model Name", accessorKey: "name" },
        { header: "Model Type", accessorKey: "model_type"},
        {
            id: "details",
            cell: ({ row }) => (
                <Link to={`/models/${row.getValue("name")}-${row.getValue("model_type")}`.toLowerCase()} className="">
                    <ArrowRight />
                </Link>
            ),
        },
    ];

    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                <h1 className="text-xl font-bold">Available models</h1>
                <p>The following models are available for classification and filtering.</p>
                <Separator className="my-4" />
                <DataTable
                    columns={columns}
                    data={tableData || []}
                />
            </div>
        </main>
    )
}