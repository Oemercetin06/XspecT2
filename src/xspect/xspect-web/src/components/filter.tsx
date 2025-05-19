import FilterForm from "./filter-form";
import { Separator } from "./ui/separator";

export default function Filter() {
    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                <h1 className="text-xl font-bold">Filter your data</h1>
                <p>Upload your data and select parameters to filter it.</p>
                <Separator className="my-4" />
                <FilterForm />
            </div>
        </main>
    )
}