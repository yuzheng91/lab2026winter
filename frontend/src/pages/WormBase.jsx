import React, { useState, useMemo } from "react";
import {
  Box,
  Button,
  TextField,
  Typography,
  Paper,
  CircularProgress,
  Alert,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Grid,
  FormGroup,
  FormControlLabel,
  Checkbox,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Chip,
  Divider,
} from "@mui/material";
import { DataGrid, GridToolbar } from "@mui/x-data-grid";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import SendIcon from "@mui/icons-material/Send";
import ManageAccountsIcon from "@mui/icons-material/ManageAccounts";
import TableChartIcon from "@mui/icons-material/TableChart";
import DescriptionIcon from "@mui/icons-material/Description";
import DownloadIcon from "@mui/icons-material/Download";
import BarChartIcon from "@mui/icons-material/BarChart"; // Graphic View Icon

// [新增] 引入 Recharts 圖表套件
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  Label,
  ZAxis,
  LabelList,
  Cell,
} from "recharts";

import { analyzeWormBase } from "../api/wormBase";

const ASPECTS = [
  { key: "F", label: "Molecular Function (MF)", btnLabel: "GO Term (MF)" },
  { key: "P", label: "Biological Process (BP)", btnLabel: "GO Term (BP)" },
  { key: "C", label: "Cellular Component (CC)", btnLabel: "GO Term (CC)" },
];

// --- 欄位定義 ---
const getColumns = (methodLabel) => [
  { field: "GO ID", headerName: "GO ID", flex: 1, minWidth: 130 },
  {
    field: "direction",
    headerName: "Type",
    width: 120,
    renderCell: (params) => {
      const isEnriched = params.value === "Enriched";
      return (
        <Chip
          label={params.value}
          size="small"
          sx={{
            bgcolor: isEnriched ? "#ffebee" : "#e3f2fd",
            color: isEnriched ? "#c62828" : "#1565c0",
            fontWeight: "bold",
            borderRadius: 1,
            height: 24,
          }}
        />
      );
    },
  },
  {
    field: "Expected",
    headerName: "Expected",
    flex: 1.2,
    minWidth: 120,
    valueGetter: (value, row) => {
      const r = row || value?.row;
      if (!r) return "-";
      const count = r["Gene Count"] || 0;
      const total = r["TotalInAspect"] || 1;
      const ratio = count / total;
      return `${count}/${total} (${(ratio * 100).toFixed(2)}%)`;
    },
  },
  {
    field: "Observed",
    headerName: "Observed",
    flex: 1.2,
    minWidth: 120,
    valueGetter: (value, row) => {
      const r = row || value?.row;
      if (!r) return "-";
      const hit = r["Input Hit"] || 0;
      const total = r["TotalInput"] || 1;
      const ratio = hit / total;
      return `${hit}/${total} (${(ratio * 100).toFixed(2)}%)`;
    },
  },
  {
    field: "Fold Change",
    headerName: "Fold Change",
    type: "number",
    flex: 0.8,
    minWidth: 110,
    valueFormatter: (value) => (value == null ? "" : Number(value).toFixed(2)),
  },
  {
    field: "display_p",
    headerName: `${methodLabel} P-value`,
    type: "number",
    flex: 1,
    minWidth: 130,
    valueFormatter: (value) =>
      value == null ? "" : Number(value).toExponential(3),
  },
];

// --- 表格組件 ---
const AspectResultPanel = ({ data, methodLabel }) => {
  if (!data || data.length === 0)
    return (
      <Typography sx={{ p: 4, textAlign: "center", color: "#666" }}>
        No significant terms found.
      </Typography>
    );

  const columns = getColumns(methodLabel);

  return (
    <Box sx={{ mt: 2, height: 600, width: "100%" }}>
      <DataGrid
        rows={data}
        columns={columns}
        initialState={{
          pagination: { paginationModel: { pageSize: 10 } },
          sorting: { sortModel: [{ field: "display_p", sort: "asc" }] },
        }}
        pageSizeOptions={[10, 25, 50, 100]}
        disableRowSelectionOnClick
        slots={{ toolbar: GridToolbar }}
        density="compact"
        getRowId={(row) => row["GO ID"]}
      />
    </Box>
  );
};

// --- 公式組件 ---
const FoldEnrichmentFormula = () => (
  <Box sx={{ display: "flex", alignItems: "center", gap: 1, my: 1 }}>
    <Typography variant="body2" fontWeight="bold" color="text.primary">
      Fold Enrichment:
    </Typography>
    <Box
      sx={{
        display: "inline-flex",
        flexDirection: "column",
        alignItems: "center",
        verticalAlign: "middle",
      }}
    >
      <Typography
        variant="caption"
        sx={{
          borderBottom: "1px solid black",
          px: 0.5,
          lineHeight: 1.2,
          fontSize: "0.85rem",
        }}
      >
        Observed Ratio
      </Typography>
      <Typography
        variant="caption"
        sx={{ px: 0.5, lineHeight: 1.2, fontSize: "0.85rem" }}
      >
        Expected Ratio
      </Typography>
    </Box>
  </Box>
);

// --- [新增] Graphic View 組件 ---
const GraphicViewPanel = ({ data, threshold }) => {
  if (!data || data.length === 0)
    return (
      <Typography sx={{ p: 4, textAlign: "center" }}>
        No data for visualization
      </Typography>
    );

  // 1. 準備 Volcano Plot 資料
  // X: Fold Change (Backend 已經是 log2)
  // Y: -log10(P-value)
  const volcanoData = data.map((row) => ({
    ...row,
    x: row["Fold Change"],
    y: -Math.log10(row.display_p || 1e-10), // 避免 log(0)
    color: row.direction === "Enriched" ? "#ef5350" : "#42a5f5", // 紅 vs 藍
  }));

  // 計算 Y 軸最大值，讓 cutoff 線漂亮一點
  const maxY = Math.max(...volcanoData.map((d) => d.y), 0) + 1;
  // 計算 Cutoff 線 (-log10 of threshold)
  const cutoffLine = -Math.log10(parseFloat(threshold) || 0.05);

  // 2. 準備 Bubble Plot 資料
  // 為了讓 Y 軸顯示 GO Term (文字)，我們需要 mapping
  // 我們取前 20 個最顯著的 (避免圖太長)
  const bubbleData = data.map((row, index) => ({
    ...row,
    x: -Math.log10(row.display_p || 1e-10), // X: log2(Fold Enrichment)
    y: index, // Y: 使用 index 來排列
    z: row["Input Hit"], // Size: Gene Count
    name: row["GO ID"], // Label
    color: row.direction === "Enriched" ? "#ef5350" : "#42a5f5",
  }));

  // Custom Tooltip
  const CustomTooltip = ({ active, payload }) => {
    if (active && payload && payload.length) {
      const d = payload[0].payload;
      return (
        <Paper sx={{ p: 1, bgcolor: "rgba(255,255,255,0.9)" }}>
          <Typography variant="body2" fontWeight="bold">
            {d["GO ID"]}
          </Typography>
          <Typography variant="caption" display="block">
            Type: {d.direction}
          </Typography>
          <Typography variant="caption" display="block">
            P-value: {d.display_p?.toExponential(2)}
          </Typography>
          <Typography variant="caption" display="block">
            log2(Fold): {d["Fold Change"]?.toFixed(2)}
          </Typography>
          <Typography variant="caption" display="block">
            Count: {d["Input Hit"]}
          </Typography>
        </Paper>
      );
    }
    return null;
  };

  return (
    <Box sx={{ p: 2 }}>
      {/* --- Volcano Plot --- */}
      <Typography
        variant="subtitle1"
        align="center"
        gutterBottom
        fontWeight="bold"
      >
        Volcano Plot
      </Typography>
      <Typography
        variant="caption"
        align="center"
        display="block"
        color="text.secondary"
        sx={{ mb: 2 }}
      >
        X: log2(Fold Enrichment), Y: -log10(P-value)
      </Typography>

      <Box sx={{ height: 350, width: "100%", mb: 4 }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis
              type="number"
              dataKey="x"
              name="log2(Fold)"
              label={{
                value: "log2(Fold Enrichment)",
                position: "bottom",
                offset: 0,
              }}
            />
            <YAxis
              type="number"
              dataKey="y"
              name="-log10(P)"
              label={{
                value: "-log10(P-value)",
                angle: -90,
                position: "insideLeft",
              }}
              domain={[0, "auto"]}
            />
            <Tooltip content={<CustomTooltip />} />
            <ReferenceLine
              y={cutoffLine}
              stroke="red"
              strokeDasharray="3 3"
              label="Cutoff"
            />
            <Scatter name="Genes" data={volcanoData} fill="#8884d8">
              {volcanoData.map((entry, index) => (
                <cell key={`cell-${index}`} fill={entry.color} />
              ))}
            </Scatter>
          </ScatterChart>
        </ResponsiveContainer>
      </Box>

      <Divider sx={{ my: 4 }} />

      {/* --- Bubble Plot --- */}
      <Typography
        variant="subtitle1"
        align="center"
        gutterBottom
        fontWeight="bold"
      >
        Bubble Plot
      </Typography>
      <Typography
        variant="caption"
        align="center"
        display="block"
        color="text.secondary"
        sx={{ mb: 2 }}
      >
        X: -log10(P-value), Size: Gene Count
      </Typography>

      <Box sx={{ height: 400, width: "100%" }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 100 }}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis
              type="number"
              dataKey="x"
              label={{
                value: "-log10(P-value)",
                position: "bottom",
                offset: 0,
              }}
            />
            {/* Y軸設為 Category，並用 tickFormatter 顯示 GO ID */}
            <YAxis
              type="number"
              dataKey="y"
              name="GO Term"
              tickCount={data.length}
              tickFormatter={(val) => {
                const item = bubbleData.find((d) => d.y === val);
                return item ? item.name : "";
              }}
              domain={[0, "auto"]}
              width={90} // 增加寬度以顯示文字
              interval={0}
            />
            <ZAxis type="number" dataKey="z" range={[100, 2000]} name="Count" />
            <Tooltip
              content={<CustomTooltip />}
              cursor={{ strokeDasharray: "3 3" }}
            />
            <Scatter name="GO Terms" data={bubbleData} fill="#8884d8">
              {bubbleData.map((entry, index) => (
                <cell key={`cell-${index}`} fill={entry.color} />
              ))}
              <LabelList
                dataKey="x"
                position="center"
                formatter={(val) => Number(val).toFixed(1)} // 取小數點 1 位 (例如 1.7)
                style={{
                  fill: "#fff", // 白色文字
                  fontSize: "8px", // 字體大小
                  fontWeight: "bold", // 粗體
                  pointerEvents: "none", // 避免擋住 Tooltip
                }}
              />
            </Scatter>
          </ScatterChart>
        </ResponsiveContainer>
      </Box>
    </Box>
  );
};

// --- Main Component ---
export default function WormBase() {
  const [listName, setListName] = useState("");
  const [inputText, setInputText] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [rawResult, setRawResult] = useState(null);
  const [method, setMethod] = useState("FDR");
  const [threshold, setThreshold] = useState("0.01");
  const [selectedAspects, setSelectedAspects] = useState({
    F: true,
    P: true,
    C: true,
  });
  const [tabValue, setTabValue] = useState("P");

  const handleAspectChange = (event) => {
    const { name, checked } = event.target;
    setSelectedAspects((prev) => {
      const newState = { ...prev, [name]: checked };
      if (tabValue === name && !checked) {
        const firstAvailable = ASPECTS.find((a) => newState[a.key])?.key;
        if (firstAvailable) setTabValue(firstAvailable);
      }
      if (!ASPECTS.find((a) => prev[a.key]) && checked) {
        setTabValue(name);
      }
      return newState;
    });
  };

  const handleAnalyze = async () => {
    if (!inputText.trim()) return;
    if (!Object.values(selectedAspects).some((v) => v)) {
      setError("Please select at least one aspect.");
      return;
    }
    setLoading(true);
    setError(null);
    setRawResult(null);
    try {
      const data = await analyzeWormBase(inputText, "None", "1.0");
      setRawResult(data);
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
    setListName("Sample_Worm_Gene_List");
    const sampleGenes = [
      "WBGene00002324",
      "WBGene00003164",
      "WBGene00004343",
      "WBGene00007012",

      "WBGene00007013",
      "WBGene00007017",
      "WBGene00007018",
      "WBGene00007026",

      "WBGene00020951",
      "WBGene00001078",
      "WBGene00007189",
      "WBGene00010720",

      "WBGene00017983",
      "WBGene00018175",
      "WBGene00020375",
      "WBGene00020820",

      "WBGene00003582",
      "WBGene00004879",
      "WBGene00004880",
      "WBGene00004881",

      "WBGene00004882",
      "WBGene00004883",
      "WBGene00004884",
      "WBGene00004885",

      "WBGene00010551",
      "WBGene00012343",
      "WBGene00012904",
      "WBGene00013188",

      "WBGene00015943",
      "WBGene00018156",
      "WBGene00021365",
      "WBGene00021929",
    ];
    setInputText(sampleGenes.join("\n"));
  };

  const handleClear = () => {
    setInputText("");
    setListName("");
  };

  const displayResult = useMemo(() => {
    if (!rawResult) return null;
    let pKey = "p_raw";
    if (method === "FDR") pKey = "p_fdr";
    if (method === "Bonferroni") pKey = "p_bon";
    const threshVal = threshold === "None" ? 1.0 : parseFloat(threshold);

    const processList = (list, direction) => {
      if (!list) return [];
      return list
        .filter((row) => row[pKey] <= threshVal)
        .map((row) => ({
          ...row,
          display_p: row[pKey],
          direction: direction,
        }));
    };

    const mergeData = (aspectData) => {
      if (!aspectData) return [];
      const enriched = processList(aspectData.enriched, "Enriched");
      const depleted = processList(aspectData.depleted, "Depleted");
      return [...enriched, ...depleted].sort(
        (a, b) => a.display_p - b.display_p,
      );
    };

    return {
      meta: rawResult.meta,
      F: mergeData(rawResult.F),
      P: mergeData(rawResult.P),
      C: mergeData(rawResult.C),
    };
  }, [rawResult, method, threshold]);

  const styles = {
    sectionBorder: "1px solid #ddd",
    headerBg: "#f9f9f9",
    accentBlueBg: "#e3f2fd",
    accentBlueText: "#1976d2",
    btnGreen: "#1b5e20",
    outerBorder: "1px solid #ccc",
    mainHeaderBg: "#dcdcdc",
    subHeaderBg: "#f5f5f5",
    tableBorderColor: "#e0e0e0",
    labelColumnBg: "#f0f4f8",
    darkBarBg: "#1c384e",
  };

  // 輔助函式：取得目前的 Data
  const getCurrentData = () => {
    if (!displayResult) return [];
    if (tabValue === "F" && selectedAspects.F) return displayResult.F;
    if (tabValue === "P" && selectedAspects.P) return displayResult.P;
    if (tabValue === "C" && selectedAspects.C) return displayResult.C;
    return [];
  };

  return (
    <Box sx={{ p: 4, maxWidth: 1400, margin: "0 auto" }}>
      <Typography variant="h5" sx={{ fontWeight: "bold", mb: 3 }}>
        YMLA: C. elegans Gene Enrichment Analysis
      </Typography>

      <Grid container spacing={3}>
        {/* Step 1 Input */}
        <Grid size={{ xs: 12, md: 5 }}>
          <Paper
            elevation={0}
            sx={{ border: styles.sectionBorder, height: "100%" }}
          >
            <Box
              sx={{
                p: 2,
                borderBottom: styles.sectionBorder,
                bgcolor: styles.headerBg,
              }}
            >
              <Typography variant="subtitle1" fontWeight="bold">
                Step1. Input gene list
              </Typography>
            </Box>
            <Box sx={{ p: 2 }}>
              <Box sx={{ display: "flex", gap: 1, mb: 2 }}>
                <TextField
                  size="small"
                  fullWidth
                  placeholder="List Name"
                  value={listName}
                  onChange={(e) => setListName(e.target.value)}
                />
                <Button
                  variant="contained"
                  size="small"
                  onClick={handleClear}
                  sx={{
                    bgcolor: "#757575",
                    color: "white",
                    "&:hover": { bgcolor: "#616161" },
                  }}
                >
                  Clear all
                </Button>
              </Box>
              <TextField
                multiline
                rows={10}
                fullWidth
                variant="outlined"
                value={inputText}
                onChange={(e) => setInputText(e.target.value)}
                sx={{ mb: 2, bgcolor: "#fff" }}
                placeholder="Paste IDs here..."
              />
              <Button
                variant="contained"
                color="primary"
                sx={{ textTransform: "none" }}
                onClick={handleSample}
              >
                Load Sample Input
              </Button>
            </Box>
          </Paper>
        </Grid>

        {/* Step 2 Features */}
        <Grid size={{ xs: 12, md: 7 }}>
          <Paper
            elevation={0}
            sx={{ border: styles.sectionBorder, height: "100%" }}
          >
            <Box
              sx={{
                p: 2,
                borderBottom: styles.sectionBorder,
                bgcolor: styles.headerBg,
              }}
            >
              <Typography variant="subtitle1" fontWeight="bold">
                Step2. Select biological feature(s)
              </Typography>
            </Box>
            <Box sx={{ p: 2 }}>
              <Accordion
                defaultExpanded
                elevation={0}
                sx={{ border: "1px solid #bbdefb", mb: 2 }}
              >
                <AccordionSummary
                  expandIcon={
                    <ExpandMoreIcon sx={{ color: styles.accentBlueText }} />
                  }
                  sx={{
                    bgcolor: styles.accentBlueBg,
                    minHeight: 48,
                    "& .MuiAccordionSummary-content": { margin: "12px 0" },
                  }}
                >
                  <Typography color={styles.accentBlueText} fontWeight="bold">
                    Ontology and Annotation Features
                  </Typography>
                </AccordionSummary>
                <AccordionDetails sx={{ pt: 2 }}>
                  <Typography
                    variant="body2"
                    sx={{
                      textDecoration: "underline",
                      color: "primary.main",
                      mb: 1,
                    }}
                  >
                    Gene Ontology (GO)
                  </Typography>
                  <FormGroup row sx={{ ml: 1 }}>
                    {ASPECTS.map((aspect) => (
                      <FormControlLabel
                        key={aspect.key}
                        control={
                          <Checkbox
                            size="small"
                            checked={selectedAspects[aspect.key]}
                            onChange={handleAspectChange}
                            name={aspect.key}
                          />
                        }
                        label={
                          <Typography variant="body2">
                            {aspect.label}
                          </Typography>
                        }
                        sx={{ width: "45%" }}
                      />
                    ))}
                  </FormGroup>
                </AccordionDetails>
              </Accordion>
              <Box sx={{ mt: 3, p: 2, bgcolor: "#f5f5f5", borderRadius: 1 }}>
                <Typography
                  variant="caption"
                  fontWeight="bold"
                  display="block"
                  mb={1}
                >
                  Analysis Parameters
                </Typography>
                <Box sx={{ display: "flex", gap: 2 }}>
                  <FormControl
                    size="small"
                    sx={{ minWidth: 150, bgcolor: "white" }}
                  >
                    <InputLabel>Correction</InputLabel>
                    <Select
                      value={method}
                      label="Correction"
                      onChange={(e) => setMethod(e.target.value)}
                    >
                      <MenuItem value="None">None</MenuItem>
                      <MenuItem value="FDR">FDR</MenuItem>
                      <MenuItem value="Bonferroni">Bonferroni</MenuItem>
                    </Select>
                  </FormControl>
                  <FormControl
                    size="small"
                    sx={{ minWidth: 150, bgcolor: "white" }}
                  >
                    <InputLabel>P-value Cutoff</InputLabel>
                    <Select
                      value={threshold}
                      label="P-value Cutoff"
                      onChange={(e) => setThreshold(e.target.value)}
                    >
                      <MenuItem value="0.05">0.05</MenuItem>
                      <MenuItem value="0.01">0.01</MenuItem>
                      <MenuItem value="None">No Cutoff</MenuItem>
                    </Select>
                  </FormControl>
                </Box>
              </Box>
            </Box>
          </Paper>
        </Grid>
      </Grid>

      <Box sx={{ mt: 4, display: "flex", justifyContent: "center" }}>
        <Button
          variant="contained"
          size="large"
          onClick={handleAnalyze}
          disabled={loading}
          startIcon={!loading && <SendIcon />}
          sx={{
            minWidth: 250,
            bgcolor: styles.btnGreen,
            fontSize: "1.1rem",
            textTransform: "none",
            py: 1.5,
            "&:hover": { bgcolor: "#144a17" },
          }}
        >
          {loading ? (
            <CircularProgress size={26} color="inherit" />
          ) : (
            "Run Enrichment Analysis"
          )}
        </Button>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mt: 2, mx: "auto", maxWidth: 800 }}>
          {error}
        </Alert>
      )}

      {displayResult && (
        <Paper
          elevation={0}
          sx={{ mt: 5, border: styles.outerBorder, overflow: "hidden" }}
        >
          {/* Main Header */}
          <Box
            sx={{
              bgcolor: styles.mainHeaderBg,
              p: 1.5,
              display: "flex",
              alignItems: "center",
              borderBottom: styles.outerBorder,
            }}
          >
            <DescriptionIcon sx={{ mr: 1 }} />
            <Typography variant="h6" fontWeight="bold">
              Enrichment Analysis Result
            </Typography>
          </Box>

          <Box sx={{ p: 2 }}>
            {/* 2. User's Specification */}
            <Paper elevation={0} sx={{ border: styles.outerBorder, mb: 3 }}>
              <Box
                sx={{
                  bgcolor: styles.subHeaderBg,
                  p: 1,
                  display: "flex",
                  alignItems: "center",
                  borderBottom: styles.tableBorderColor,
                }}
              >
                <ManageAccountsIcon sx={{ mr: 1 }} />
                <Typography variant="subtitle1" fontWeight="bold">
                  User's Specification
                </Typography>
              </Box>
              <Box>
                <Box
                  sx={{
                    display: "flex",
                    borderBottom: `1px solid ${styles.tableBorderColor}`,
                  }}
                >
                  <Box
                    sx={{
                      flex: "0 0 30%",
                      bgcolor: styles.labelColumnBg,
                      p: 1,
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "flex-end",
                      borderRight: `1px solid ${styles.tableBorderColor}`,
                    }}
                  >
                    <Typography variant="body2" fontWeight="bold">
                      Name of each input list
                    </Typography>
                  </Box>
                  <Box sx={{ flex: 1, p: 1, pl: 2 }}>
                    <Typography variant="body2">
                      {listName || "Translation Efficiency Mid"}
                    </Typography>
                  </Box>
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    borderBottom: `1px solid ${styles.tableBorderColor}`,
                  }}
                >
                  <Box
                    sx={{
                      flex: "0 0 30%",
                      bgcolor: styles.labelColumnBg,
                      p: 1,
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "flex-end",
                      borderRight: `1px solid ${styles.tableBorderColor}`,
                    }}
                  >
                    <Typography variant="body2" fontWeight="bold">
                      # of genes in each input list
                    </Typography>
                  </Box>
                  <Box sx={{ flex: 1, p: 1, pl: 2 }}>
                    <Typography variant="body2">
                      {displayResult.meta.mapped_live}
                    </Typography>
                  </Box>
                </Box>
                <Box
                  sx={{
                    bgcolor: styles.darkBarBg,
                    p: 0.5,
                    textAlign: "center",
                  }}
                >
                  <Typography variant="caption" color="white">
                    See the analysis result of a chosen feature
                  </Typography>
                </Box>
                <Box
                  sx={{
                    display: "flex",
                    borderTop: `1px solid ${styles.tableBorderColor}`,
                  }}
                >
                  <Box
                    sx={{
                      flex: "0 0 30%",
                      bgcolor: styles.labelColumnBg,
                      p: 2,
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "flex-end",
                      borderRight: `1px solid ${styles.tableBorderColor}`,
                    }}
                  >
                    <Typography variant="body2" fontWeight="bold">
                      Ontology and Annotation Features
                    </Typography>
                  </Box>
                  <Box sx={{ flex: 1, p: 2, display: "flex", gap: 1 }}>
                    {ASPECTS.map((aspect) => {
                      const isActive = tabValue === aspect.key;
                      if (!selectedAspects[aspect.key]) return null;
                      return (
                        <Button
                          key={aspect.key}
                          variant="contained"
                          onClick={() => setTabValue(aspect.key)}
                          disableElevation
                          sx={{
                            textTransform: "none",
                            boxShadow: "none",
                            fontSize: "0.85rem",
                            px: 2,
                            border: "1px solid",
                            borderColor: isActive ? "#0d6efd" : "#999",
                            bgcolor: isActive ? "#0d6efd" : "#e0e0e0",
                            color: isActive ? "#fff" : "#000",
                            "&:hover": {
                              bgcolor: isActive ? "#0b5ed7" : "#d5d5d5",
                              borderColor: isActive ? "#0b5ed7" : "#888",
                            },
                          }}
                        >
                          {aspect.btnLabel}
                        </Button>
                      );
                    })}
                  </Box>
                </Box>
              </Box>
            </Paper>

            {/* 3. Table View */}
            <Paper elevation={0} sx={{ border: styles.outerBorder, mb: 3 }}>
              <Box
                sx={{
                  bgcolor: styles.subHeaderBg,
                  p: 1,
                  display: "flex",
                  alignItems: "center",
                  borderBottom: styles.tableBorderColor,
                }}
              >
                <TableChartIcon sx={{ mr: 1 }} />
                <Typography variant="subtitle1" fontWeight="bold">
                  Table View
                </Typography>
              </Box>
              <Box sx={{ p: 2 }}>
                <Box sx={{ mb: 2, p: 2, bgcolor: "#fafafa", borderRadius: 1 }}>
                  <Typography variant="body2">
                    <strong>A:</strong> # of input genes with this GO term
                  </Typography>
                  <Typography variant="body2">
                    <strong>B:</strong> # of input genes
                  </Typography>
                  <Typography variant="body2">
                    <strong>C:</strong> # of genes (in the <em>C. elegans</em>{" "}
                    genome) with this GO term
                  </Typography>
                  <Typography variant="body2">
                    <strong>D:</strong> # of genes in the <em>C. elegans</em>{" "}
                    genome ({displayResult.meta.total_background} genes)
                  </Typography>
                  <FoldEnrichmentFormula />
                </Box>
                <Box sx={{ mb: 2 }}>
                  <Button
                    variant="outlined"
                    startIcon={<DownloadIcon />}
                    size="small"
                    sx={{
                      textTransform: "none",
                      color: "#333",
                      borderColor: "#ccc",
                    }}
                  >
                    Download as a .csv file
                  </Button>
                </Box>
                <Box sx={{ minHeight: 400 }}>
                  {!Object.values(selectedAspects).some((v) => v) ? (
                    <Typography sx={{ p: 4, textAlign: "center" }}>
                      Please select an aspect above.
                    </Typography>
                  ) : (
                    <>
                      {tabValue === "F" && selectedAspects.F && (
                        <AspectResultPanel
                          data={displayResult.F}
                          methodLabel={method === "None" ? "Raw" : method}
                        />
                      )}
                      {tabValue === "P" && selectedAspects.P && (
                        <AspectResultPanel
                          data={displayResult.P}
                          methodLabel={method === "None" ? "Raw" : method}
                        />
                      )}
                      {tabValue === "C" && selectedAspects.C && (
                        <AspectResultPanel
                          data={displayResult.C}
                          methodLabel={method === "None" ? "Raw" : method}
                        />
                      )}
                    </>
                  )}
                </Box>
              </Box>
            </Paper>

            {/* 4. Graphic View (新加入區塊) */}
            <Paper elevation={0} sx={{ border: styles.outerBorder }}>
              <Box
                sx={{
                  bgcolor: styles.subHeaderBg,
                  p: 1,
                  display: "flex",
                  alignItems: "center",
                  borderBottom: styles.tableBorderColor,
                }}
              >
                <BarChartIcon sx={{ mr: 1 }} />
                <Typography variant="subtitle1" fontWeight="bold">
                  Graphic View
                </Typography>
              </Box>
              <GraphicViewPanel data={getCurrentData()} threshold={threshold} />
            </Paper>
          </Box>
        </Paper>
      )}
    </Box>
  );
}
