// src/api/wormBase.js
const API_BASE_URL = 'http://127.0.0.1:8000'; // 請依實際狀況調整

export const analyzeWormBase = async (text, method, threshold) => {
  const response = await fetch(`${API_BASE_URL}/worm/api/analyze/`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ text, method, threshold }),
  });

  if (!response.ok) {
    const errorData = await response.json();
    throw new Error(errorData.error || 'Network response was not ok');
  }

  return response.json();
};