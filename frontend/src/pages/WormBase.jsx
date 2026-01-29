import React, { useState, useEffect } from "react";
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
  Stack,
  Grid,
  FormGroup,
  FormControlLabel,
  Checkbox,
} from "@mui/material";
import { DataGrid, GridToolbar } from "@mui/x-data-grid";
import { analyzeWormBase } from "../api/wormBase";
import BugReportIcon from "@mui/icons-material/BugReport";

// 共用的欄位定義
const columns = [
  { field: "GO ID", headerName: "GO ID", flex: 1, minWidth: 150 },
  { field: "Expected", headerName: "Expected Ratio", flex: 1.2, minWidth: 180 },
  { field: "Observed", headerName: "Observed Ratio", flex: 1.2, minWidth: 180 },
  {
    field: "Fold Change",
    headerName: "Fold Change",
    type: "number",
    flex: 0.8,
    minWidth: 120,
    valueFormatter: (value) => (value == null ? "" : Number(value).toFixed(2)),
  },
  {
    field: "P-value",
    headerName: "Adj. P-value",
    type: "number",
    flex: 1,
    minWidth: 130,
    valueFormatter: (value) =>
      value == null ? "" : Number(value).toExponential(3),
  },
];

// 定義 Aspect 的對照表，方便程式碼迴圈生成
const ASPECTS = [
  { key: "F", label: "Molecular Function (F)" },
  { key: "P", label: "Biological Process (P)" },
  { key: "C", label: "Cellular Component (C)" },
];

const AspectResultPanel = ({ data, aspectName }) => {
  if (!data) return null;

  return (
    <Box sx={{ mt: 2 }}>
      <Typography variant="h6" color="primary" gutterBottom>
        {aspectName} - Enriched
      </Typography>
      <Box sx={{ height: 400, width: "100%", mb: 4 }}>
        <DataGrid
          rows={data.enriched.map((r, i) => ({ ...r, id: r["GO ID"] || i }))}
          columns={columns}
          initialState={{ pagination: { paginationModel: { pageSize: 5 } } }}
          pageSizeOptions={[5, 10, 25]}
          disableRowSelectionOnClick
          slots={{ toolbar: GridToolbar }}
        />
      </Box>

      <Typography variant="h6" color="secondary" gutterBottom>
        {aspectName} - Depleted
      </Typography>
      <Box sx={{ height: 400, width: "100%" }}>
        <DataGrid
          rows={data.depleted.map((r, i) => ({ ...r, id: r["GO ID"] || i }))}
          columns={columns}
          initialState={{ pagination: { paginationModel: { pageSize: 5 } } }}
          pageSizeOptions={[5, 10, 25]}
          disableRowSelectionOnClick
          slots={{ toolbar: GridToolbar }}
        />
      </Box>
    </Box>
  );
};

export default function WormBase() {
  const [inputText, setInputText] = useState("");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);

  // 參數 State
  const [method, setMethod] = useState("FDR");
  const [threshold, setThreshold] = useState("0.05");

  // 【新增】Aspect 勾選狀態 (預設全選)
  const [selectedAspects, setSelectedAspects] = useState({
    F: true,
    P: true,
    C: true,
  });

  // Tab 狀態 (改用 key 字串 'F', 'P', 'C' 來控制，比較安全)
  const [tabValue, setTabValue] = useState("P");

  // 處理 Aspect 勾選變化
  const handleAspectChange = (event) => {
    const { name, checked } = event.target;
    setSelectedAspects((prev) => {
      const newState = { ...prev, [name]: checked };

      // 如果目前顯示的 Tab 被取消勾選了，自動切換到第一個有被勾選的 Tab
      if (tabValue === name && !checked) {
        const firstAvailable = ASPECTS.find((a) => newState[a.key])?.key;
        if (firstAvailable) setTabValue(firstAvailable);
      }
      // 如果目前沒有選中任何 Tab (例如之前全取消)，且現在勾選了一個，就切過去
      if (!ASPECTS.find((a) => prev[a.key]) && checked) {
        setTabValue(name);
      }
      return newState;
    });
  };

  const handleAnalyze = async () => {
    if (!inputText.trim()) return;

    // 防呆：如果使用者全部取消勾選
    if (!Object.values(selectedAspects).some((v) => v)) {
      setError("Please select at least one aspect (F, P, or C) to analyze.");
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);

    try {
      // 這裡目前還是呼叫後端跑全部 (後端通常很快)，前端只負責隱藏顯示
      // 如果資料量真的很大，也可以把 selectedAspects 傳給後端只跑部分
      const data = await analyzeWormBase(inputText, method, threshold);
      setResult(data);

      // 分析完成後，確保 Tab 停留在第一個有勾選且有資料的分頁
      const firstAvailable = ASPECTS.find((a) => selectedAspects[a.key])?.key;
      if (firstAvailable) setTabValue(firstAvailable);
    } catch (err) {
      console.error(err);
      setError(err.message || "Analysis failed.");
    } finally {
      setLoading(false);
    }
  };

  const handleSample = () => {
    setInputText(
      "abl-1\natl-1\nclk-2\nhim-6\nhpr-9\nhpr-17\nhsr-9\nhus-1\nmrt-2\nD1053.2\ndot-1.2\ndot-1.4\ndot-1.3\nuri-1\ngen-1\ndot-1.1\ndot-1.5"
    );
  };

  return (
    <Box sx={{ p: 4, maxWidth: 1200, margin: "0 auto" }}>
      <Box sx={{ display: "flex", alignItems: "center", mb: 1 }}>
        <Typography variant="h5" sx={{ fontWeight: "bold" }}>
          YMLA: Yeast Multiple List Analyzer
        </Typography>
      </Box>
      <Box sx={{ borderBottom: "1px solid #333", mb: 3 }} />

      {/* 輸入與設定區塊 (改用 Grid 做左右佈局) */}
      <Paper
        elevation={3}
        sx={{ p: 3, mb: 4, bgcolor: "white", border: "1px solid #bbb" }}
      >
        <Grid container spacing={3}>
          {/* 左邊：輸入框 (佔 8/12) */}
          <Grid item xs={12} md={8}>
            <Box
              sx={{ display: "flex", justifyContent: "space-between", mb: 1 }}
            >
              <Typography variant="subtitle1" fontWeight="bold">
                Gene List Input
              </Typography>
              <Button size="small" onClick={handleSample}>
                Load Sample
              </Button>
            </Box>
            <TextField
              multiline
              rows={8}
              fullWidth
              variant="outlined"
              value={inputText}
              onChange={(e) => setInputText(e.target.value)}
              placeholder="Paste C. elegans Gene Names or WBGeneIDs here..."
              sx={{ bgcolor: "#fafafa" }}
            />
          </Grid>

          {/* 右邊：參數設定 (佔 4/12) */}
          <Grid item xs={12} md={4}>
            <Box
              sx={{
                display: "flex",
                flexDirection: "column",
                height: "100%",
                gap: 2,
              }}
            >
              <Typography variant="subtitle1" fontWeight="bold">
                Parameters
              </Typography>

              {/* 參數選擇 */}
              <FormControl size="small" fullWidth>
                <InputLabel>Correction method</InputLabel>
                <Select
                  value={method}
                  label="Correction method"
                  onChange={(e) => setMethod(e.target.value)}
                >
                  <MenuItem value="None">None</MenuItem>
                  <MenuItem value="FDR">FDR</MenuItem>
                  <MenuItem value="Bonferroni">Bonferroni</MenuItem>
                </Select>
              </FormControl>

              <FormControl size="small" fullWidth>
                <InputLabel>Corrected p-value cutoff</InputLabel>
                <Select
                  value={threshold}
                  label="Corrected p-value cutoff"
                  onChange={(e) => setThreshold(e.target.value)}
                >
                  <MenuItem value="0.05">0.05</MenuItem>
                  <MenuItem value="0.01">0.01</MenuItem>
                  <MenuItem value="0.001">0.001</MenuItem>
                  <MenuItem value="0.000001">1e-6</MenuItem>
                </Select>
              </FormControl>

              {/* Aspect 勾選區 */}
              <Box
                sx={{ border: "1px solid #ddd", borderRadius: 1, p: 2, mt: 1 }}
              >
                <Typography variant="caption" color="text.secondary">
                  Select Aspects to Display:
                </Typography>
                <FormGroup>
                  <FormControlLabel
                    control={
                      <Checkbox
                        checked={selectedAspects.F}
                        onChange={handleAspectChange}
                        name="F"
                      />
                    }
                    label="Molecular Function (F)"
                  />
                  <FormControlLabel
                    control={
                      <Checkbox
                        checked={selectedAspects.P}
                        onChange={handleAspectChange}
                        name="P"
                      />
                    }
                    label="Biological Process (P)"
                  />
                  <FormControlLabel
                    control={
                      <Checkbox
                        checked={selectedAspects.C}
                        onChange={handleAspectChange}
                        name="C"
                      />
                    }
                    label="Cellular Component (C)"
                  />
                </FormGroup>
              </Box>

              {/* 按鈕置底 */}
            </Box>
          </Grid>
        </Grid>

        {error && (
          <Alert severity="error" sx={{ mt: 2 }}>
            {error}
          </Alert>
        )}
      </Paper>
      <Box sx={{ my: 4, display: "flex", justifyContent: "center" }}>
        <Button
          variant="contained"
          size="large"
          onClick={handleAnalyze}
          disabled={loading}
          sx={{ minWidth: 200}}
        >
          {loading ? (
            <CircularProgress size={24} color="inherit" />
          ) : (
            "Run Analysis"
          )}
        </Button>
      </Box>
      {/* 結果區 */}
      {result && (
        <Paper sx={{ width: "100%", p: 2 }}>
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Mapped {result.meta.mapped_live} live genes out of{" "}
              {result.meta.input_total} inputs.
            </Typography>
          </Box>

          {/* 上方 Tabs：只顯示被勾選的 Aspect */}
          {/* 注意：這裡 Tab 的 value 我們直接用 'F', 'P', 'C' 字串 */}
          <Tabs
            value={tabValue}
            onChange={(e, v) => setTabValue(v)}
            variant="scrollable"
            scrollButtons="auto"
            indicatorColor="primary"
            textColor="primary"
            sx={{ borderBottom: 1, borderColor: "divider" }}
          >
            {/* 根據勾選狀態動態產生 Tab */}
            {selectedAspects.F && (
              <Tab label="Molecular Function (F)" value="F" />
            )}
            {selectedAspects.P && (
              <Tab label="Biological Process (P)" value="P" />
            )}
            {selectedAspects.C && (
              <Tab label="Cellular Component (C)" value="C" />
            )}
          </Tabs>

          {/* 內容區：根據 Tab 顯示對應的 AspectResultPanel */}
          <Box sx={{ p: 2, bgcolor: "#f9f9f9", minHeight: 500 }}>
            {!Object.values(selectedAspects).some((v) => v) ? (
              <Typography
                sx={{ p: 4, textAlign: "center", color: "text.secondary" }}
              >
                Please select at least one aspect to view results.
              </Typography>
            ) : (
              <>
                {tabValue === "F" && selectedAspects.F && (
                  <AspectResultPanel
                    data={result.F}
                    aspectName="Molecular Function"
                  />
                )}
                {tabValue === "P" && selectedAspects.P && (
                  <AspectResultPanel
                    data={result.P}
                    aspectName="Biological Process"
                  />
                )}
                {tabValue === "C" && selectedAspects.C && (
                  <AspectResultPanel
                    data={result.C}
                    aspectName="Cellular Component"
                  />
                )}
              </>
            )}
          </Box>
        </Paper>
      )}
    </Box>
  );
}
