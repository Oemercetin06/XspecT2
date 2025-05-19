import axios from 'axios';
import { Models, ClassificationResult, ModelMetadata } from '@/types';

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