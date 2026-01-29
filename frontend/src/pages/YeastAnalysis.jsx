import React, { useState, useEffect } from 'react';
import { 
  Box, 
  Button, 
  TextField, 
  Typography, 
  Paper, 
  Tab, 
  Tabs,
  CircularProgress,
  Alert,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Link,
  Stack
} from '@mui/material';
import { DataGrid, GridToolbar } from '@mui/x-data-grid';
import SearchIcon from '@mui/icons-material/Search'; // 記得安裝 icons: npm install @mui/icons-material
import DeleteIcon from '@mui/icons-material/Delete';
import { analyzeYeast } from '../api/yeastApi';

// --- DataGrid 欄位定義 (保持不變) ---
const columns = [
  { field: 'Domain Name', headerName: 'Domain Name', flex: 1.5, minWidth: 200 },
  { field: 'Expected Ratio', headerName: 'Expected Ratio', flex: 1.2, minWidth: 180 },
  { field: 'Observed Ratio', headerName: 'Observed Ratio', flex: 1.2, minWidth: 180 },
  { 
    field: 'Fold Change', headerName: 'Fold Change', type: 'number', flex: 0.8, minWidth: 120,
    valueFormatter: (value) => value == null ? '' : Number(value).toFixed(2)
  },
  { 
    field: 'P-value', headerName: 'P-value', type: 'number', flex: 1, minWidth: 130,
    valueFormatter: (value) => value == null ? '' : Number(value).toExponential(3)
  },
  { 
    field: 'FDR', headerName: 'FDR', type: 'number', flex: 1, minWidth: 130,
    valueFormatter: (value) => value == null ? '' : Number(value).toExponential(3)
  },
  { 
    field: 'Bonferroni', headerName: 'Bonferroni', type: 'number', flex: 1, minWidth: 130,
    valueFormatter: (value) => value == null ? '' : Number(value).toExponential(3)
  },
];

// 範例資料
const SAMPLE_GENES = "YAL001C\nYAL002W\nYAL003W\nYBL004W\nYCR005C";

export default function YeastAnalysis() {
  // --- State ---
  const [inputText, setInputText] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [tabValue, setTabValue] = useState(0);
  
  const [rawResult, setRawResult] = useState(null);
  const [filteredResult, setFilteredResult] = useState(null);

  // 篩選條件
  const [method, setMethod] = useState('None');
  const [threshold, setThreshold] = useState('0.05');

  // --- Functions ---
  const handleAnalyze = async () => {
    if (!inputText.trim()) {
        setError("Please enter at least one gene.");
        return;
    }
    
    setLoading(true);
    setError(null);
    try {
      const data = await analyzeYeast(inputText, "None", "None");
      
      const enrichedWithId = data.enriched.map((row) => ({ ...row, id: row['Domain Name'] }));
      const depletedWithId = data.depleted.map((row) => ({ ...row, id: row['Domain Name'] }));
      
      setRawResult({ enriched: enrichedWithId, depleted: depletedWithId, meta: data.meta });
      setTabValue(0);
    } catch (err) {
      console.error(err);
      setError("Analysis failed. Please check your input or server status.");
    } finally {
      setLoading(false);
    }
  };

  const handleReset = () => {
    setInputText('');
    setRawResult(null);
    setFilteredResult(null);
    setError(null);
  };

  // 填入範例功能
  const fillSample = () => {
    setInputText(SAMPLE_GENES);
    setError(null);
  };

  // 前端過濾邏輯
  useEffect(() => {
    if (!rawResult) return;

    const doFilter = (dataList) => {
      if (!threshold || threshold === 'None') return dataList;
      const threshVal = parseFloat(threshold);
      
      let colName = 'P-value';
      if (method === 'FDR') colName = 'FDR';
      if (method === 'Bonferroni') colName = 'Bonferroni';

      return dataList.filter(row => {
        const val = row[colName];
        return val !== null && val !== undefined && val <= threshVal;
      });
    };

    const newEnriched = doFilter(rawResult.enriched);
    const newDepleted = doFilter(rawResult.depleted);

    setFilteredResult({
      enriched: newEnriched,
      depleted: newDepleted,
      meta: {
        ...rawResult.meta,
        enriched_count: newEnriched.length,
        depleted_count: newDepleted.length
      }
    });

  }, [rawResult, method, threshold]);

  return (
    <Box sx={{ p: 4, maxWidth: 1200, margin: '0 auto' }}>
      
      {/* 標題區域 */}
      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
        <SearchIcon sx={{ fontSize: 32, mr: 1 }} />
        <Typography variant="h5" sx={{ fontWeight: 'bold', color: 'black' }}>
          YPIBP: A Repository for Yeast PI-Binding Proteins
        </Typography>
      </Box>

      <Box sx={{ borderBottom: '1px solid #333', mb: 3 }} />

      {/* 輸入區塊容器 (灰色背景) */}
      <Paper 
        elevation={3} 
        sx={{ 
            p: 3, 
            mb: 4, 
            bgcolor: '#e0e0e0', // 模仿目標圖的灰色背景
            border: '1px solid #bbb'
        }}
      >
        {/* Sample Input 連結區 */}
        <Box sx={{ mb: 1 }}>
            <Typography variant="subtitle1" component="span" sx={{ fontWeight: 'bold', mr: 1 }}>
                Sample Input:
            </Typography>
            <Box component="div">
                <Link 
                    component="button" 
                    variant="body2" 
                    onClick={fillSample}
                    sx={{ textAlign: 'left', display: 'block', mb: 0.5 }}
                >
                    Input a list of proteins
                </Link>
            </Box>
        </Box>

        {/* 純白色的輸入框 */}
        <TextField
          multiline
          rows={8}
          fullWidth
          variant="outlined"
          value={inputText}
          onChange={(e) => setInputText(e.target.value)}
          placeholder="Paste Gene List here..."
          sx={{ 
            bgcolor: 'white', // 輸入框保持白色
            '& .MuiOutlinedInput-root': {
                padding: 1, // 調整內距讓文字不要貼邊
            }
          }}
        />

        {/* 按鈕區域 (底部置中) */}
        <Stack direction="row" spacing={2} justifyContent="center" sx={{ mt: 3 }}>
            <Button 
                variant="contained" 
                color="primary" // 藍色 Search
                size="large" 
                onClick={handleAnalyze}
                disabled={loading}
                sx={{ minWidth: 120, fontWeight: 'bold' }}
            >
                {loading ? <CircularProgress size={24} color="inherit" /> : "Search"}
            </Button>
            
            <Button 
                variant="contained" 
                color="warning" // 黃色 Reset
                size="large" 
                onClick={handleReset}
                sx={{ minWidth: 120, fontWeight: 'bold', color: 'white' }}
            >
                Reset
            </Button>
        </Stack>

        {error && <Alert severity="error" sx={{ mt: 2 }}>{error}</Alert>}
      </Paper>

      {/* 篩選與結果顯示區 (維持原本設計，但稍微美化) */}
      {filteredResult && (
        <Box>
            {/* 篩選工具列 */}
            <Paper variant="outlined" sx={{ p: 2, mb: 2, bgcolor: '#f5f5f5', display: 'flex', gap: 2, alignItems: 'center' }}>
                <Typography variant="subtitle2" sx={{ fontWeight: 'bold' }}>Filter Results:</Typography>
                <FormControl size="small" sx={{ minWidth: 150, bgcolor: 'white' }}>
                    <InputLabel>Method</InputLabel>
                    <Select value={method} label="Method" onChange={(e) => setMethod(e.target.value)}>
                        <MenuItem value="None">None</MenuItem>
                        <MenuItem value="FDR">FDR</MenuItem>
                        <MenuItem value="Bonferroni">Bonferroni</MenuItem>
                    </Select>
                </FormControl>
                <FormControl size="small" sx={{ minWidth: 150, bgcolor: 'white' }}>
                    <InputLabel>Threshold</InputLabel>
                    <Select value={threshold} label="Threshold" onChange={(e) => setThreshold(e.target.value)}>
                        <MenuItem value="0.05">0.05</MenuItem>
                        <MenuItem value="0.01">0.01</MenuItem>
                        <MenuItem value="0.001">0.001</MenuItem>
                        <MenuItem value="None">None</MenuItem>
                    </Select>
                </FormControl>
                
                <Typography variant="body2" sx={{ ml: 'auto', color: 'text.secondary' }}>
                   Found: <b>{filteredResult.meta.enriched_count}</b> enriched domains
                </Typography>
            </Paper>

            {/* 結果表格 */}
            <Paper sx={{ width: '100%', borderTop: '3px solid #1976d2' }}> {/* 頂部加一條藍線增加質感 */}
            <Tabs 
                value={tabValue} 
                onChange={(e, v) => setTabValue(v)} 
                indicatorColor="primary"
                textColor="primary"
                sx={{ borderBottom: 1, borderColor: 'divider' }}
            >
                <Tab label={`Enriched Domains (${filteredResult.enriched.length})`} />
                <Tab label={`Depleted Domains (${filteredResult.depleted.length})`} />
            </Tabs>

            <Box sx={{ height: 600, width: '100%', p: 0 }}>
                <DataGrid
                rows={tabValue === 0 ? filteredResult.enriched : filteredResult.depleted}
                columns={columns}
                initialState={{
                    pagination: { paginationModel: { pageSize: 10 } },
                    sorting: { sortModel: [{ field: 'P-value', sort: 'asc' }] },
                }}
                pageSizeOptions={[10, 25, 50, 100]}
                disableRowSelectionOnClick
                slots={{ toolbar: GridToolbar }}
                slotProps={{ toolbar: { showQuickFilter: true } }}
                />
            </Box>
            </Paper>
        </Box>
      )}
    </Box>
  );
}