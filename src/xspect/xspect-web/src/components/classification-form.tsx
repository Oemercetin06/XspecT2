import { zodResolver } from "@hookform/resolvers/zod"
import { useForm } from "react-hook-form"
import { z } from "zod"
import { Button } from "@/components/ui/button"
import { Switch } from "@/components/ui/switch"
import { Input } from "@/components/ui/input"
import {
    Form,
    FormControl,
    FormDescription,
    FormField,
    FormItem,
    FormLabel,
    FormMessage,
} from "@/components/ui/form"
import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from "@/components/ui/select"
import { classify, getModels, uploadFile } from "../api"
import {
    FileUpload,
    FileUploadDropzone,
    FileUploadItem,
    FileUploadItemDelete,
    FileUploadItemMetadata,
    FileUploadItemPreview,
    FileUploadItemProgress,
    FileUploadList,
    FileUploadTrigger,
} from "@/components/ui/file-upload";
import { Upload, X } from "lucide-react";
import { useState, useCallback, useEffect } from "react"
import { useNavigate } from "react-router-dom"

const FormSchema = z.object({
    input_file: z.string().min(1, "Please upload a file"),
    classification_type: z.enum([
        "Genus",
        "Species",
        "MLST",
    ]),
    model: z.string(),
    sparse_sampling: z.boolean(),
    sparse_sampling_step: z.number().min(1).max(500)
});

export default function ClassificationForm() {

    const [models, setModels] = useState<Record<string, string[]>>({});

    useEffect(() => {
        getModels()
            .then((data) => {
                setModels(data);
            })
            .catch((error) => {
                console.error("Error fetching models:", error);
            });
    }, []);

    const navigate = useNavigate()

    const form = useForm<z.infer<typeof FormSchema>>({
        resolver: zodResolver(FormSchema),
        defaultValues: {
            classification_type: "Species",
            model: undefined,
            sparse_sampling: false,
            sparse_sampling_step: 1,
        },
    })
    const [files, setFiles] = useState<File[]>([]);


    const onUpload = useCallback(
        async (
            files: File[],
            {
                onSuccess,
                onError,
            }: {
                onSuccess: (file: File) => void;
                onError: (file: File, error: Error) => void;
            },
        ) => {
            const file = files[0]
            try {
                const { filename } = await uploadFile(file)
                onSuccess(file)
                console.log("File uploaded successfully:", filename)
                form.setValue("input_file", filename)
            } catch (error) {
                onError(file, error as Error)
                console.error("Error uploading file:", error)
            }
        },
        [form]
    );

    const onFileReject = useCallback((file: File, message: string) => {
        console.log(message, {
            description: `"${file.name.length > 20 ? `${file.name.slice(0, 20)}...` : file.name}" has been rejected`,
        });
    }, []);

    function onSubmit(data: z.infer<typeof FormSchema>) {
        classify(data.input_file, data.classification_type, data.model, data.sparse_sampling_step)
            .then((response) => {
                console.log("Classification response:", response)
                navigate(`/result/${response.uuid}`)
            })
            .catch((error) => {
                console.error("Error during classification:", error)
            });
    }

    return (

        <Form {...form}>
            <form onSubmit={form.handleSubmit(onSubmit)} className="space-y-6">
                <FormField
                    control={form.control}
                    name="input_file"
                    render={({ field }) => (
                        <FormItem>
                            <FormLabel>Input File</FormLabel>
                            <FormControl>
                                <FileUpload
                                    value={files}
                                    onValueChange={setFiles}
                                    onUpload={onUpload}
                                    onFileReject={onFileReject}
                                    maxFiles={1}
                                    multiple
                                >
                                    {files.length === 0 && (
                                        <FileUploadDropzone>
                                            <div className="flex flex-col items-center gap-1 text-center">
                                                <div className="flex items-center justify-center rounded-full border p-2.5">
                                                    <Upload className="size-6 text-muted-foreground" />
                                                </div>
                                                <p className="font-medium text-sm">Drag & drop file here</p>
                                                <p className="text-muted-foreground text-xs">
                                                    Or click to browse
                                                </p>
                                            </div>
                                            <FileUploadTrigger asChild>
                                                <Button variant="outline" size="sm" className="mt-2 w-fit">
                                                    Browse files
                                                </Button>
                                            </FileUploadTrigger>
                                        </FileUploadDropzone>
                                    )}
                                    <FileUploadList>
                                        {files.map((file, index) => (
                                            <FileUploadItem key={index} value={file} className="flex-col">
                                                <div className="flex w-full items-center gap-2">
                                                    <FileUploadItemPreview />
                                                    <FileUploadItemMetadata />
                                                    <FileUploadItemDelete asChild>
                                                        <Button variant="ghost" size="icon" className="size-7" onClick={() => {
                                                            setFiles((prev) => prev.filter((_, i) => i !== index));
                                                            form.setValue("input_file", "")
                                                        }}>
                                                            <X />
                                                        </Button>
                                                    </FileUploadItemDelete>
                                                </div>
                                                <FileUploadItemProgress />
                                            </FileUploadItem>
                                        ))}
                                    </FileUploadList>
                                </FileUpload>
                            </FormControl>
                            <FormDescription>
                                Upload the file you would like to classify.
                            </FormDescription>
                            <FormMessage />
                        </FormItem>
                    )}
                />
                <FormField
                    control={form.control}
                    name="classification_type"
                    render={({ field }) => (
                        <FormItem>
                            <FormLabel>Classification Type</FormLabel>
                            <FormControl>
                                <Select defaultValue="Species" onValueChange={(value) => field.onChange(value)}>
                                    <SelectTrigger className="w-full">
                                        <SelectValue />
                                    </SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="Genus">Genus</SelectItem>
                                        <SelectItem value="Species">Species</SelectItem>
                                        <SelectItem value="MLST">MLST</SelectItem>
                                    </SelectContent>
                                </Select>
                            </FormControl>
                            <FormDescription>
                                Select the type of classification you would like to perform.
                            </FormDescription>
                            <FormMessage />
                        </FormItem>
                    )}
                />
                <FormField
                    control={form.control}
                    name="model"
                    render={({ field }) => (
                        <FormItem>
                            <FormLabel>Model</FormLabel>
                            <FormControl>
                                <Select onValueChange={field.onChange}>
                                    <SelectTrigger className="w-full">
                                        <SelectValue />
                                    </SelectTrigger>
                                    <SelectContent>
                                        {models[form.watch("classification_type")]?.map((model: string) => (
                                            <SelectItem key={model} value={model}>
                                                {model}
                                            </SelectItem>
                                        ))}
                                    </SelectContent>
                                </Select>
                            </FormControl>
                            <FormDescription>
                                Select the model you would like to use for classification.
                            </FormDescription>
                            <FormMessage />
                        </FormItem>
                    )}
                />
                <FormField
                    control={form.control}
                    name="sparse_sampling"
                    render={({ field }) => (
                        <FormItem className="flex flex-row items-center justify-between pr-2">
                            <div className="space-y-0.5">
                                <FormLabel>Sparse Sampling</FormLabel>
                                <FormDescription>
                                    Enable sparse sampling for classification.
                                </FormDescription>
                            </div>
                            <FormControl>
                                <Switch
                                    checked={field.value}
                                    onCheckedChange={field.onChange}
                                />
                            </FormControl>
                        </FormItem>
                    )}
                />
                {form.watch("sparse_sampling") && (
                    <FormField
                        control={form.control}
                        name="sparse_sampling_step"
                        render={({ field }) => (
                            <FormItem>
                                <FormLabel>Sparse Sampling Step</FormLabel>
                                <FormControl>
                                    <Input
                                        type="number"
                                        value={field.value}
                                        onChange={(e) => field.onChange(Number(e.target.value))}
                                    />
                                </FormControl>
                                <FormDescription>
                                    Set the step size for sparse sampling.
                                </FormDescription>
                                <FormMessage />
                            </FormItem>
                        )}
                    />
                )}
                <Button type="submit">Classify</Button>
            </form>
        </Form>
    )
}
