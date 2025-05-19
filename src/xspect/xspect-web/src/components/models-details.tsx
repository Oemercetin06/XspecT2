import { ModelMetadata } from "@/types";
import { useEffect, useState } from "react";
import { getModelMetadata } from "../api";
import { Link, useParams } from "react-router-dom";
import { Separator } from "@/components/ui/separator";
import { DataTable } from "@/components/data-table"

export default function ModelDetails() {
    const { model_slug } = useParams();
    const [modelMetadata, setModelMetadata] = useState<ModelMetadata | null>(null);

    useEffect(() => {
        if (model_slug) {
            getModelMetadata(model_slug)
                .then((data: ModelMetadata) => {
                    setModelMetadata(data);
                })
                .catch((err) => console.error('Fetch error:', err));
        }
    }, [model_slug]);

    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                {modelMetadata && (
                    <>
                        <h1 className="text-xl font-bold">
                            {modelMetadata.model_display_name} {modelMetadata.model_type} Model
                        </h1>

                        <div className="flex h-5 items-center space-x-4">
                            <div>k = {modelMetadata.k}</div>
                            <Separator orientation="vertical" />
                            <div>num_hashes = {modelMetadata.num_hashes}</div>
                            <Separator orientation="vertical" />
                            <div>fpr = {modelMetadata.fpr}</div>

                            {modelMetadata.model_type === "Species" && (
                                <>
                                    <Separator orientation="vertical" />
                                    <div>kernel = {modelMetadata.kernel}</div>
                                    <Separator orientation="vertical" />
                                    <div>C = {modelMetadata.C}</div>
                                </>
                            )}
                        </div>

                        <Separator className="my-4" />

                        <h2 className="text-l font-bold">Author</h2>
                        <div className="flex flex-col gap-2">
                            <p>Model author: {modelMetadata.author == null ? "N/A" : modelMetadata.author}</p>
                            <p>
                                Model author email: {modelMetadata.author_email == null
                                    ? "N/A"
                                    : <Link to={`mailto:${modelMetadata.author_email}`} className="font-medium underline underline-offset-4">
                                        {modelMetadata.author_email}
                                    </Link>
                                }
                            </p>
                        </div>

                        <Separator className="my-4" />

                        <h2 className="text-l font-bold mb-2">Model display names</h2>
                        <DataTable
                            columns={[
                                { accessorKey: "key", header: "ID" },
                                { accessorKey: "value", header: "Display name" }
                            ]}
                            data={Object.entries(modelMetadata.display_names).map(([key, value]) => ({ key, value }))}
                        />

                        <Separator className="my-4" />

                        <h2 className="text-l font-bold mb-2">K-mer Index Training Data</h2>
                        {typeof modelMetadata.training_accessions === "object" && Array.isArray(modelMetadata.training_accessions) && (
                            <p>
                                The model was trained on the following accessions:
                                <ul className="list-disc list-inside">
                                    {modelMetadata.training_accessions.map((accession, i) => (
                                        <li key={i}>
                                            <Link
                                                to={`https://www.ncbi.nlm.nih.gov/datasets/genome/${accession}/`}
                                                target="_blank"
                                                rel="noopener noreferrer"
                                                className="font-medium underline underline-offset-4"
                                            >
                                                {accession}
                                            </Link>
                                        </li>
                                    ))}
                                </ul>
                            </p>
                        )}
                        {typeof modelMetadata.training_accessions === "object" && !Array.isArray(modelMetadata.training_accessions) && (
                            <DataTable
                                columns={[
                                    { accessorKey: "key", header: "Name" },
                                    {
                                        accessorKey: "value", header: "Accessions", cell: ({ getValue }) => {
                                            const value = getValue();
                                            return Array.isArray(value)
                                                ? value.map((accession, i) => (
                                                    <>
                                                        {i > 0 && ", "}
                                                        <Link
                                                            to={`https://www.ncbi.nlm.nih.gov/datasets/genome/${accession}/`}
                                                            target="_blank"
                                                            rel="noopener noreferrer"
                                                            className="font-medium underline underline-offset-4"
                                                        >
                                                            {accession}
                                                        </Link>
                                                    </>
                                                ))
                                                : value;
                                        }
                                    }
                                ]}
                                data={Object.entries(modelMetadata.training_accessions).map(([key, value]) => ({
                                    key: modelMetadata.display_names[key],
                                    value
                                }))}
                            />
                        )}
                    </>
                )}
                {!modelMetadata && (
                    <>
                        <h1 className="text-xl font-bold">Error</h1>
                        <p>No model could be found with slug <span className="font-mono">{model_slug}</span>.</p>
                    </>
                )}
            </div>
        </main>
    )
}