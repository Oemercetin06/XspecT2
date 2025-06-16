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
import { filterSequences, getModelMetadata, getModels, uploadFile } from "../api"
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
import { Check, ChevronsUpDown, Upload, X } from "lucide-react";
import { useState, useCallback, useEffect } from "react"
import { useNavigate } from "react-router-dom"
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { ModelMetadata } from "../types"
import { Popover, PopoverContent, PopoverTrigger } from "./ui/popover"
import { Command, CommandEmpty, CommandGroup, CommandInput, CommandItem, CommandList } from "./ui/command"
import { cn } from "@/lib/utils"

const FormSchema = z.object({
    input_file: z.string().min(1, "Please upload a file"),
    filter_type: z.enum([
        "Genus",
        "Species",
    ]),
    model: z.string(),
    filter_species: z.string().optional(),
    filter_mode: z.enum([
        "threshold",
        "maximum"
    ]),
    threshold: z.number().min(0).max(1),
    sparse_sampling: z.boolean(),
    sparse_sampling_step: z.number().min(1).max(500)
});

export default function FilterForm() {

    const [models, setModels] = useState<Record<string, string[]>>({});
    const [modelMetadata, setModelMetadata] = useState<ModelMetadata | null>(null);

    const form = useForm<z.infer<typeof FormSchema>>({
        resolver: zodResolver(FormSchema),
        defaultValues: {
            filter_type: "Species",
            model: undefined,
            filter_mode: "threshold",
            threshold: 0.7,
            sparse_sampling: false,
            sparse_sampling_step: 1,
        },
    })
    const [files, setFiles] = useState<File[]>([]);

    useEffect(() => {
        getModels()
            .then((data) => {
                setModels(data);
            })
            .catch((error) => {
                console.error("Error fetching models:", error);
            });
    }, []);

    const selectedModel = form.watch("model") !== undefined && form.watch("filter_type") !== undefined
        ? `${form.watch("model")}-${form.watch("filter_type")}`.toLowerCase()
        : undefined;

    useEffect(() => {
        if (models && selectedModel) {
            getModelMetadata(selectedModel)
                .then((data) => {
                    setModelMetadata(data);
                })
                .catch((error) => {
                    console.error("Error fetching model metadata:", error);
                });
        }
    }, [models, selectedModel, form]);

    const navigate = useNavigate()


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
        const threshold = data.filter_mode === "threshold" ? data.threshold : -1;
        filterSequences(
            data.filter_type,
            data.model,
            data.input_file,
            threshold,
            data.filter_species,
            data.sparse_sampling_step
        ).then((response) => {
            console.log("Filtering response:", response)
            navigate(`/filter-result/${response.uuid}`)
        })
        .catch((error) => {
            console.error("Error during filtering:", error)
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
                                Upload the file you would like to filter.
                            </FormDescription>
                            <FormMessage />
                        </FormItem>
                    )}
                />
                <FormField
                    control={form.control}
                    name="filter_type"
                    render={({ field }) => (
                        <FormItem>
                            <FormLabel>Filter Type</FormLabel>
                            <FormControl>
                                <Select defaultValue="Species" onValueChange={(value) => field.onChange(value)}>
                                    <SelectTrigger className="w-full">
                                        <SelectValue />
                                    </SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="Genus">Genus</SelectItem>
                                        <SelectItem value="Species">Species</SelectItem>
                                    </SelectContent>
                                </Select>
                            </FormControl>
                            <FormDescription>
                                Select the type of filter you would like to perform.
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
                                        {models[form.watch("filter_type")]?.map((model: string) => (
                                            <SelectItem key={model} value={model}>
                                                {model}
                                            </SelectItem>
                                        ))}
                                    </SelectContent>
                                </Select>
                            </FormControl>
                            <FormDescription>
                                Select the model you would like to use for filtering.
                            </FormDescription>
                            <FormMessage />
                        </FormItem>
                    )}
                />
                {form.watch("filter_type") === "Species" && (
                    <>
                        {form.watch("model") != undefined && (
                            <FormField
                                control={form.control}
                                name="filter_species"
                                render={({ field }) => (
                                    <FormItem className="flex flex-col">
                                        <FormLabel>Filter Species</FormLabel>
                                        <Popover>
                                            <PopoverTrigger asChild>
                                                <FormControl>
                                                    <Button
                                                        variant="outline"
                                                        role="combobox"
                                                        className={cn(
                                                            "w-full justify-between",
                                                            !field.value && "text-muted-foreground"
                                                        )}
                                                    >
                                                        {field.value
                                                            ? modelMetadata?.display_names[field.value] as string
                                                            : "Select species..."}
                                                        <ChevronsUpDown className="ml-2 h-4 w-4 shrink-0 opacity-50" />
                                                    </Button>
                                                </FormControl>
                                            </PopoverTrigger>
                                            <PopoverContent className="p-0">
                                                <Command>
                                                    <CommandInput placeholder="Search species..." />
                                                    <CommandList>
                                                        <CommandEmpty>No species found.</CommandEmpty>
                                                        <CommandGroup>
                                                            {modelMetadata?.display_names && Object.entries(modelMetadata.display_names).map(([value, label]) => (
                                                                <CommandItem
                                                                    value={label as string}
                                                                    key={value}
                                                                    onSelect={() => {
                                                                        form.setValue("filter_species", value)
                                                                    }}
                                                                >
                                                                    {label as string}
                                                                    <Check
                                                                        className={cn(
                                                                            "ml-auto",
                                                                            value === field.value
                                                                                ? "opacity-100"
                                                                                : "opacity-0"
                                                                        )}
                                                                    />
                                                                </CommandItem>
                                                            ))}
                                                        </CommandGroup>
                                                    </CommandList>
                                                </Command>
                                            </PopoverContent>
                                        </Popover>
                                        <FormDescription>
                                            Select the species you would like to filter for.
                                        </FormDescription>
                                        <FormMessage />
                                    </FormItem>
                                )}
                            />
                        )}

                        <FormField
                            control={form.control}
                            name="filter_mode"
                            render={({ field }) => (
                                <FormItem>
                                    <FormLabel>Filter Mode</FormLabel>
                                    <FormControl>
                                        <Tabs defaultValue="threshold" onValueChange={field.onChange}>
                                            <TabsList className="w-full">
                                                <TabsTrigger value="threshold">Filter by threshold</TabsTrigger>
                                                <TabsTrigger value="maximum">Filter by maximum scoring species</TabsTrigger>
                                            </TabsList>
                                        </Tabs>
                                    </FormControl>
                                    <FormDescription>
                                        Select the filtering mode. You can either filter by threshold or by maximum scoring species, meaning that contigs will only be included if the selected species is the highest-scoring one.
                                    </FormDescription>
                                    <FormMessage />
                                </FormItem>
                            )}
                        />
                    </>
                )}
                {form.watch("filter_mode") === "threshold" && (
                    <FormField
                        control={form.control}
                        name="threshold"
                        render={({ field }) => (
                            <FormItem>
                                <FormLabel>Threshold</FormLabel>
                                <FormControl>
                                    <Input
                                        type="number"
                                        value={field.value}
                                        defaultValue={0.7}
                                        step={0.01}
                                        min={0}
                                        max={1}
                                        onChange={(e) => field.onChange(Number(e.target.value))}
                                    />
                                </FormControl>
                                <FormDescription>
                                    Set the threshold for filtering.
                                </FormDescription>
                                <FormMessage />
                            </FormItem>
                        )}
                    />
                )}
                <FormField
                    control={form.control}
                    name="sparse_sampling"
                    render={({ field }) => (
                        <FormItem className="flex flex-row items-center justify-between pr-2">
                            <div className="space-y-0.5">
                                <FormLabel>Sparse Sampling</FormLabel>
                                <FormDescription>
                                    Enable sparse sampling for filtering.
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
                {
                    form.watch("sparse_sampling") && (
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
                    )
                }
                <Button type="submit">Filter</Button>
            </form >
        </Form >
    )
}
