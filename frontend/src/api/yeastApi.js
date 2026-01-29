// src/api/yeastApi.js
import axios from 'axios';

// 假設你的 Django 後端跑在 8000 port
const API_BASE_URL = 'http://localhost:8000';

export const analyzeYeast = async (text) => {
  try {
    const response = await axios.post(`${API_BASE_URL}/yeast/api/analyze/`, {
      text: text
    });
    return response.data;
  } catch (error) {
    console.error("Error analyzing yeast data:", error);
    throw error;
  }
};
