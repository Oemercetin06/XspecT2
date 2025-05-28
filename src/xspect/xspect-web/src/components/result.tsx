import { Link, useParams } from "react-router-dom"
import { getClassificationResult, getModelMetadata } from "../api";
import { ModelMetadata, ClassificationResult } from "../types";
import { useState, useEffect, use } from "react";
import { LoadingSpinner } from "./spinner";
import { ResultChart, ResultChartProps } from "./result-chart";
import { Separator } from "@/components/ui/separator";
import { Button } from "./ui/button";
import { DropdownMenuCheckboxes } from "./dropdown-checkboxes";
import { DropdownMenuSlider } from "./dropdown-slider";

type CheckboxItem = {
    id: string;
    label: string;
    checked: boolean;
};

const useCheckboxItems = () => {
    const [items, setItems] = useState<CheckboxItem[]>([]);

    const handleCheckedChange = (id: string, checked: boolean) => {
        setItems(items.map(item =>
            item.id === id ? { ...item, checked } : item
        ));
    };

    const checkboxItems = items.map(item => ({
        ...item,
        onCheckedChange: (checked: boolean) => handleCheckedChange(item.id, checked)
    }));

    return { checkboxItems, setItems };
};


export default function ResultPage() {
    const { checkboxItems: contig_checkbox_items, setItems: setCheckboxItems } = useCheckboxItems();
    const { classification_uuid } = useParams();
    const [classificationResult, setClassificationResult] = useState<ClassificationResult | null>(null);
    const [chartData, setChartData] = useState<ResultChartProps[] | null>(null);
    const [modelMetadata, setModelMetadata] = useState<ModelMetadata | null>(null);
    const [numResults, setNumResults] = useState(15);

    const checkedItemsState = JSON.stringify(
        contig_checkbox_items.map(item => ({ id: item.id, checked: item.checked }))
    );

    useEffect(() => {
        const fetchResult = () => {
            if (classification_uuid) {
                getClassificationResult(classification_uuid).then((data) => {
                    setClassificationResult(data);
                    getModelMetadata(data.model_slug).then((modelMetadata) => {
                        setModelMetadata(modelMetadata);
                    }).catch((error) => {
                        console.error('Error fetching model metadata:', error);
                        setTimeout(fetchResult, 500);
                    });
                }).catch((error) => {
                    console.error('Error fetching classification result:', error);
                    setTimeout(fetchResult, 500);
                });
            }
        };
        fetchResult();
    }, [classification_uuid, setClassificationResult, setModelMetadata]);

    useEffect(() => {
        if (classificationResult) {
            const contigNames = Object.keys(classificationResult.scores);
            if (contigNames.length <= 20) {
                const initialCheckboxItems = contigNames
                    .filter(name => name !== "total")
                    .map((name) => ({
                        id: name,
                        label: name,
                        checked: true,
                    }));
                setCheckboxItems(initialCheckboxItems);
            }
        }
    }, [classificationResult, setCheckboxItems]);

    useEffect(() => {
        if (classificationResult && modelMetadata) {
            const selectedContigs = contig_checkbox_items
                .filter(item => item.checked)
                .map(item => item.id);
            console.log('selectedContigs', selectedContigs);
            const filtered_hits = Object.entries(classificationResult.hits).filter(([key]) => selectedContigs.includes(key));
            const total_hits = {}
            filtered_hits.forEach(([, value]) => {
                Object.entries(value).forEach(([label, hits]) => {
                    if (!total_hits[label]) {
                        total_hits[label] = 0;
                    }
                    total_hits[label] += hits;
                });
            });
            const total_kmers = Object.entries(classificationResult.num_kmers)
                .filter(([key]) => selectedContigs.includes(key))
                .reduce((sum, [_, value]) => sum + value, 0);
            const scores = Object.entries(total_hits).map(([label, hits]) => {
                return {
                    taxon: modelMetadata.display_names[label].replace(modelMetadata.model_display_name, modelMetadata.model_display_name.charAt(0) + '.') || label,
                    score: hits / total_kmers,
                };
            });
            scores.sort((a, b) => b.score - a.score);
            const topScores = scores.slice(0, numResults);
            setChartData(topScores);
        }
    }, [classificationResult, modelMetadata, checkedItemsState, numResults]);


    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                {!chartData && (
                    <div className="flex items-center justify-center">
                        <LoadingSpinner className="size-10" />
                        <span className="m-2 text-xl">Classifying...</span>
                    </div>
                )}
                {chartData && (<>
                    <h1 className="text-xl font-bold">Classification Results</h1>
                    <p>{classificationResult?.input_source} was classified using the <Link className="font-medium underline underline-offset-4" to={`/models/${modelMetadata?.model_slug}`}>{modelMetadata?.model_display_name} {modelMetadata?.model_type.toLowerCase()}</Link> model. {classificationResult?.sparse_sampling_step === 1 ? 'No sparse sampling was used.' : 'Sparse sampling was used with step ' + classificationResult?.sparse_sampling_step + '.'}</p>
                    <p className="mt-4"><b>Prediction:</b> {modelMetadata?.display_names[classificationResult?.prediction] || 'No prediction available.'}</p>
                    <Button className="w-full mt-4" onClick={() => {
                        const blob = new Blob([JSON.stringify(classificationResult, null, 2)], { type: 'application/json' });
                        const href = URL.createObjectURL(blob);

                        const link = document.createElement("a");
                        link.href = href;
                        link.download = "result-" + classification_uuid + ".json";
                        document.body.appendChild(link);
                        link.click();

                        document.body.removeChild(link);
                        URL.revokeObjectURL(href);
                    }}>
                        Download Full Result
                    </Button>
                    <Separator className="my-4 h-px" />
                    <div className="flex justify-end space-x-2">
                        <DropdownMenuSlider triggerButtonText="#Results" value={numResults} onValueChange={setNumResults} max={modelMetadata?.display_names.length} step={1} disabled={false} />
                        <DropdownMenuCheckboxes triggerButtonText="Contigs" labelText="Select Contigs" items={contig_checkbox_items} disabled={contig_checkbox_items.length === 0} />
                    </div>
                    <ResultChart data={chartData} />
                </>)}

            </div>
        </main>
    )
}
