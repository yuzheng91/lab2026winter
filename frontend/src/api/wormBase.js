// src/api/wormBase.js
const API_BASE_URL = 'http://127.0.0.1:8000'; // 確定是這個 port 沒錯吧！

export const analyzeWormBase = async (text, runFpc, runGeneGroup) => {
  // 加上 /worm 是對的！因為你的 config/urls.py 有規定
  const response = await fetch(`${API_BASE_URL}/worm/api/analyze/`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    // 把開關丟給後端的雙管分流器
    body: JSON.stringify({ 
      text: text, 
      run_fpc: runFpc, 
      run_gene_group: runGeneGroup 
    }),
  });

  if (!response.ok) {
    const errorData = await response.json().catch(() => ({}));
    throw new Error(errorData.error || `連線失敗，HTTP 狀態碼: ${response.status}`);
  }

  return response.json();
};

export const fetchTermEvidence = async ({ term_id, genes }) => {
  // ★ 修正這裡：加上 /worm 前綴，並統一使用 API_BASE_URL
  const response = await fetch(`${API_BASE_URL}/worm/api/fetch_term_evidence/`, { 
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ term_id, genes }),
  });

  if (!response.ok) {
    throw new Error('Failed to fetch evidence data');
  }
  return response.json();
};