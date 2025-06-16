import { Link, useParams } from "react-router-dom"
import { FilteringResult } from "../types";
import { useState, useEffect, use } from "react";
import { LoadingSpinner } from "./spinner";
import { Separator } from "@/components/ui/separator";
import { Button } from "./ui/button";
import { getFilteringResult } from "../api";
import { TriangleAlert } from "lucide-react";

export default function FilteringResultPage() {
    const { filter_uuid } = useParams();
    const [filteringResult, setFilteringResult] = useState<FilteringResult | null>(null);


    useEffect(() => {
        const fetchResult = () => {
            if (filter_uuid) {
                getFilteringResult(filter_uuid).then((data) => {
                    setFilteringResult(data);
                }).catch((error) => {
                    console.error('Error fetching filtering result:', error);
                    setTimeout(fetchResult, 500);
                });
            }
        };
        fetchResult();
    }, [filter_uuid, setFilteringResult]);


    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                {!filteringResult && (
                    <div className="flex items-center justify-center">
                        <LoadingSpinner className="size-10" />
                        <span className="m-2 text-xl">Filtering...</span>
                    </div>
                )}
                {filteringResult && (
                    <>
                        <h1 className="text-xl font-bold">Filtering Results</h1>
                        <p>Your file was filtered successfully.</p>
                        {filteringResult.message != "Filtering completed successfully." && (
                            <div className="flex items-center gap-2 text-red-500">
                                <TriangleAlert className="w-5 h-5" />
                                <p>All sequences were filtered out.</p>
                            </div>
                        )}
                        <Separator className="my-4 h-px" />
                        <div className="flex flex-col gap-3 mt-4">
                            {filteringResult.message == "Filtering completed successfully." && (
                                <Button asChild>
                                    <a href={`${window.location.origin}/api/download-filtered?uuid=${filter_uuid}`} download>Download the filtered file</a>
                                </Button>
                            )}
                            <Button variant="outline" asChild>
                                <Link to={`/result/${filter_uuid}`}>View underlying classification results</Link>
                            </Button>
                        </div>
                    </>
                )}
            </div>
        </main>
    )
}
