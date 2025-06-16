import axios from 'axios';
import { Models, ClassificationResult, FilteringResult, ModelMetadata } from '@/types';

export async function getModels(): Promise<Models> {
    try {
        const response = await axios.get('/api/list-models');
        return response.data as Models;
    } catch (error) {
        console.error('Error fetching models:', error);
        throw new Error('Failed to fetch models');
    }
}

export async function uploadFile(file: File): Promise<{ filename: string }> {
    try {
        const formData = new FormData();
        formData.append('file', file);
        const response = await axios.post('/api/upload-file', formData, {
            headers: {
                'Content-Type': 'multipart/form-data',
            },
        });
        return response.data as { filename: string };
    } catch (error) {
        console.error('Error uploading file:', error);
        throw new Error('Failed to upload file');
    }
}

export async function classify(filename: string, classification_type: string, model: string, step: number): Promise<any> {
    try {
        const response = await axios.post(`/api/classify?classification_type=${classification_type}&model=${model}&file=${filename}&step=${step}`);
        return response.data;
    } catch (error) {
        console.error('Error classifying file:', error);
        throw new Error('Failed to classify file');
    }
}

export async function getClassificationResult(classification_uuid: string): Promise<ClassificationResult> {
    try {
        const response = await axios.get(`/api/classification-result`, {
            params: { uuid: classification_uuid },
        });
        return response.data as ClassificationResult;
    } catch (error) {
        console.error('Error fetching classification result:', error);
        throw new Error('Failed to fetch classification result');
    }
}

export async function getModelMetadata(modelSlug: string): Promise<ModelMetadata> {
    try {
        const response = await axios.get(`/api/model-metadata`, {
            params: { model_slug: modelSlug },
        });
        return response.data as ModelMetadata;
    } catch (error) {
        console.error('Error fetching model metadata:', error);
        throw new Error('Failed to fetch model metadata');
    }
}

export async function filterSequences(filter_type: string, genus : string, input_file : string, threshold : number, filter_species : string = "", sparse_sampling_step : number = 1): Promise<any> {
    try {
        const response = await axios.post(`/api/filter?filter_type=${filter_type}&genus=${genus}&input_file=${input_file}&threshold=${threshold}&filter_species=${filter_species}&step=${sparse_sampling_step}`);
        return response.data;
    } catch (error) {
        console.error('Error filtering sequences:', error);
        throw new Error('Failed to filter sequences');
    }
}

export async function getFilteringResult(filter_uuid: string): Promise<FilteringResult> {
    try {
        const response = await axios.get(`/api/filtering-result`, {
            params: { uuid: filter_uuid },
        });
        return response.data as FilteringResult;
    } catch (error) {
        console.error('Error fetching filtering result:', error);
        throw new Error('Failed to fetch filtering result');
    }
}