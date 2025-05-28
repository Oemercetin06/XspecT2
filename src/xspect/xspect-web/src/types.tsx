export interface Models {
  [key: string]: string[];
}

export interface ModelTableEntry {
  "model_type": "Species" | "Genus" | "Family" | "Order" | "Class" | "Phylum" | "Kingdom";
  "name": string;
}

export interface ModelMetadata {
  model_slug: string;
  model_display_name: string;
  model_type: string;
  author: string;
  author_email: string;
  k: number;
  display_names: { [key: string]: string };
  fpr: number;
  num_hashes?: number;
  kernel?: Partial<string>;
  C?: Partial<number>;
  training_accessions: { [key: string]: string[] } | string[];
  svm_accessions: { [key: string]: string[] };
}

export interface ClassificationResult {
  input_source: string;
  model_slug: string;
  sparse_sampling_step: number;
  hits: { [subsequence: string]: { [label: string]: number } };
  scores: { [subsequence: string]: { [label: string]: number } };
  num_kmers: { [subsequence: string]: number };
  prediction: string;
}

export interface FilteringResult {
  message: string;
  uuid : string;
}