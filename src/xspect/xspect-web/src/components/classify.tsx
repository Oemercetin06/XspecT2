import ClassificationForm from "./classification-form";
import { Separator } from "@/components/ui/separator"



export default function Classify() {

    return (
        <main className="flex-1 flex flex-col items-center justify-center p-4">
            <div className="w-1/2">
                <h1 className="text-xl font-bold">Classify your data</h1>
                <p>Upload your data and select a model to classify it.</p>
                <Separator className="my-4" />
                <ClassificationForm />
            </div>
        </main>
    )
}
