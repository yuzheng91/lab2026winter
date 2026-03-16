// src/api/wqfcApi.js
import axios from 'axios';

// 假設你的 Django 後端跑在 8000 port
const API_BASE_URL = 'http://localhost:8000';

// 執行主要的比對分析 (POST)
export const runWqfcPipeline = async (l1Text, l2Text, features) => {
  try {
    const response = await axios.post(`${API_BASE_URL}/wqfc/run_pipeline/`, {
      L1: l1Text,
      L2: l2Text,
      features: features
    });
    return response.data;
  } catch (error) {
    console.error("Error running WQFC pipeline:", error);
    throw error;
  }
};

// 獲取特定 Feature 的基因詳細清單 (GET)
export const fetchWqfcGeneDetail = async (listType, feature) => {
  try {
    const response = await axios.get(`${API_BASE_URL}/wqfc/gene_detail/`, {
      params: {
        list: listType,
        feature: feature
      }
    });
    return response.data;
  } catch (error) {
    console.error("Error fetching WQFC gene details:", error);
    throw error;
  }
};