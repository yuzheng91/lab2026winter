import React, { useState, useMemo } from "react";
import {
  Box,
  Button,
  TextField,
  Typography,
  Paper,
  CircularProgress,
  Alert,
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
  Dialog,
  DialogTitle,
  DialogContent,
  IconButton,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
} from "@mui/material";
import { DataGrid, GridToolbar } from "@mui/x-data-grid";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import SendIcon from "@mui/icons-material/Send";
import ManageAccountsIcon from "@mui/icons-material/ManageAccounts";
import DownloadIcon from "@mui/icons-material/Download";
import CheckIcon from "@mui/icons-material/Check";
import SettingsIcon from "@mui/icons-material/Settings";
import RefreshIcon from "@mui/icons-material/Refresh";
import BarChartIcon from "@mui/icons-material/BarChart";
import CloseIcon from "@mui/icons-material/Close";
import GroupIcon from "@mui/icons-material/Group";
import PushPinIcon from "@mui/icons-material/PushPin";

// Recharts
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  ZAxis,
  LabelList,
  Cell,
} from "recharts";

// ★ 確保你的 API 檔案中有匯出這兩個函式
import { analyzeWormBase, fetchTermEvidence } from "../api/wormBase";

// --- 定義分析類別 ---
const ASPECTS = [
  { key: "F", label: "Molecular Function (MF)", btnLabel: "GO Term (MF)" },
  { key: "P", label: "Biological Process (BP)", btnLabel: "GO Term (BP)" },
  { key: "C", label: "Cellular Component (CC)", btnLabel: "GO Term (CC)" },
];
const FGG_FEATURE = {
  key: "G",
  label: "Gene Groups (55 Features)",
  btnLabel: "Gene Groups",
};

// --- 欄位定義 ---
const getColumns = (methodLabel, onObservedClick) => [
  {
    field: "Term Name",
    headerName: "Feature Name",
    flex: 1,
    minWidth: 250,
    renderCell: (params) => {
      const id = params.row["Term ID"];
      const name = params.row["Term Name"];
      if (id && id.startsWith("GO:") && name && name !== id) {
        return `${id} (${name})`;
      }
      return name || id;
    },
  },
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
      if (r["Expected Ratio Text"]) return r["Expected Ratio Text"];
      const count = r["Gene Count"] || 0;
      const total = r["TotalInAspect"] || 1;
      return `${count}/${total} (${((count / total) * 100).toFixed(3)}%)`;
    },
  },
  {
    field: "Observed",
    headerName: "Observed",
    flex: 1.2,
    minWidth: 120,
    renderCell: (params) => {
      const r = params.row;
      if (!r) return "-";

      const fullText =
        r["Observed Ratio Text"] ||
        `${r["Input Hit"] || 0}/${r["TotalInput"] || 1} (${(((r["Input Hit"] || 0) / (r["TotalInput"] || 1)) * 100).toFixed(3)}%)`;

      const parts = fullText.split("/");
      const numerator = parts[0];
      const restOfText = parts.length > 1 ? `/${parts.slice(1).join("/")}` : "";

      return (
        <Typography variant="body2" component="div">
          <Box
            component="span"
            sx={{
              color: "primary.main",
              cursor: "pointer",
              textDecoration: "underline",
              fontWeight: "bold",
            }}
            onClick={(e) => {
              e.stopPropagation();
              if (onObservedClick) onObservedClick(r);
            }}
          >
            {numerator}
          </Box>
          <Box component="span" sx={{ color: "text.primary" }}>
            {restOfText}
          </Box>
        </Typography>
      );
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

// --- Evidence 欄位定義 ---
const evidenceColumns = [
  { field: "WormBase_ID", headerName: "WormBase ID", flex: 1, minWidth: 150 },
  { field: "Symbol", headerName: "Symbol", flex: 1, minWidth: 120 },
  { field: "GO_ID", headerName: "GO ID", flex: 1, minWidth: 120 },
  { field: "Reference", headerName: "Reference", flex: 1.5, minWidth: 180 },
  { field: "Evidence_Code", headerName: "Evidence Code", flex: 1, minWidth: 120 },
];

// --- 表格組件 ---
const AspectResultPanel = ({ data, methodLabel, onObservedClick }) => {
  if (!data || data.length === 0)
    return (
      <Typography sx={{ p: 4, textAlign: "center", color: "#666" }}>
        No significant terms found.
      </Typography>
    );

  return (
    <Box sx={{ mt: 2, height: 600, width: "100%" }}>
      <DataGrid
        rows={data}
        columns={getColumns(methodLabel, onObservedClick)}
        initialState={{
          pagination: { paginationModel: { pageSize: 10 } },
          sorting: { sortModel: [{ field: "display_p", sort: "asc" }] },
        }}
        pageSizeOptions={[10, 25, 50, 100]}
        disableRowSelectionOnClick
        slots={{ toolbar: GridToolbar }}
        density="compact"
        getRowId={(row) => row["Term ID"]}
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
    <Box sx={{ display: "inline-flex", flexDirection: "column", alignItems: "center", verticalAlign: "middle" }}>
      <Typography variant="caption" sx={{ borderBottom: "1px solid black", px: 0.5, lineHeight: 1.2, fontSize: "0.85rem" }}>
        Observed Ratio
      </Typography>
      <Typography variant="caption" sx={{ px: 0.5, lineHeight: 1.2, fontSize: "0.85rem" }}>
        Expected Ratio
      </Typography>
    </Box>
  </Box>
);

// --- Graphic View 組件 ---
const GraphicViewPanel = ({ data, threshold }) => {
  if (!data || data.length === 0) return <Typography sx={{ p: 4, textAlign: "center" }}>No data for visualization</Typography>;

  const volcanoData = data.map((row) => ({
    ...row,
    x: row["log2_fc"] > -100 ? row["log2_fc"] : 0,
    y: row["neg_log10_p"],
    color: row.direction === "Enriched" ? "#ef5350" : "#42a5f5",
  }));

  const cutoffLine = threshold === "None" ? 0 : -Math.log10(parseFloat(threshold) || 0.05);

  const bubbleData = data.map((row, index) => ({
    ...row,
    x: row["log2_fc"] > -100 ? row["log2_fc"] : 0,
    y: index,
    z: Number(row["Fold Change"]) || 0,
    name: row["Term ID"]?.startsWith("GO:") ? `${row["Term ID"]} (${row["Term Name"]})` : row["Term Name"],
    color: row.direction === "Enriched" ? "#ef5350" : "#42a5f5",
  }));

  const CustomTooltip = ({ active, payload }) => {
    if (active && payload && payload.length) {
      const d = payload[0].payload;
      const displayName = d["Term ID"]?.startsWith("GO:") ? `${d["Term ID"]} (${d["Term Name"]})` : d["Term Name"];
      return (
        <Paper sx={{ p: 1, bgcolor: "rgba(255,255,255,0.9)", maxWidth: 300 }}>
          <Typography variant="body2" fontWeight="bold">{displayName}</Typography>
          <Typography variant="caption" display="block">Type: {d.direction}</Typography>
          <Typography variant="caption" display="block">P-value: {d.display_p?.toExponential(2)}</Typography>
          <Typography variant="caption" display="block">log2(Fold): {d["log2_fc"]?.toFixed(2)}</Typography>
          {d["Input Hit"] && <Typography variant="caption" display="block">Count: {d["Input Hit"]}</Typography>}
        </Paper>
      );
    }
    return null;
  };

  return (
    <Box sx={{ p: 2 }}>
      <Typography variant="subtitle1" align="center" gutterBottom fontWeight="bold">Volcano Plot</Typography>
      <Typography variant="caption" align="center" display="block" color="text.secondary" sx={{ mb: 2 }}>
        X: log2(Fold Enrichment), Y: -log10(P-value)
      </Typography>
      <Box sx={{ height: 350, width: "100%", mb: 4 }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis type="number" dataKey="x" name="log2(Fold)" label={{ value: "log2(Fold Enrichment)", position: "bottom", offset: 0 }} />
            <YAxis type="number" dataKey="y" name="-log10(P)" width={160} label={{ value: "-log10(P-value)", angle: -90, position: "insideLeft" }} domain={[0, "auto"]} />
            <Tooltip content={<CustomTooltip />} />
            {threshold !== "None" && <ReferenceLine y={cutoffLine} stroke="red" strokeDasharray="3 3" label="Cutoff" />}
            <Scatter name="Genes" data={volcanoData}>
              {volcanoData.map((entry, index) => <Cell key={`cell-${index}`} fill={entry.color} />)}
            </Scatter>
          </ScatterChart>
        </ResponsiveContainer>
      </Box>

      <Divider sx={{ my: 4 }} />

      <Typography variant="subtitle1" align="center" gutterBottom fontWeight="bold">Bubble Plot</Typography>
      <Typography variant="caption" align="center" display="block" color="text.secondary" sx={{ mb: 2 }}>
        X: log2(Fold Enrichment), Size: Fold Change
      </Typography>
      <Box sx={{ height: 400, width: "100%" }}>
        <ResponsiveContainer>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 0 }}>
            <CartesianGrid strokeDasharray="3 3" />
            <XAxis type="number" dataKey="x" label={{ value: "log2(Fold Enrichment)", position: "bottom", offset: 0 }} />
            <YAxis type="number" dataKey="y" tickCount={data.length} domain={[0, "auto"]} width={170} interval={0} tick={{ fontSize: 11 }}
              tickFormatter={(val) => {
                const item = bubbleData.find((d) => d.y === val);
                return item ? (item.name.length > 25 ? item.name.substring(0, 25) + "..." : item.name) : "";
              }}
            />
            <ZAxis type="number" dataKey="z" range={[100, 2000]} name="Fold Change" />
            <Tooltip content={<CustomTooltip />} cursor={{ strokeDasharray: "3 3" }} />
            <Scatter name="Features" data={bubbleData}>
              {bubbleData.map((entry, index) => <Cell key={`cell-${index}`} fill={entry.color} />)}
              <LabelList dataKey="x" position="center" formatter={(val) => Number(val).toFixed(1)} style={{ fill: "#fff", fontSize: "10px", fontWeight: "bold", pointerEvents: "none" }} />
            </Scatter>
          </ScatterChart>
        </ResponsiveContainer>
      </Box>
    </Box>
  );
};

// --- Main Component ---
export default function WormBase() {
  const [listName, setListName] = useState("Sample_Gene_List");
  const [inputText, setInputText] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const [rawResult, setRawResult] = useState(null);
  const [selectedAspects, setSelectedAspects] = useState({ F: true, P: true, C: true, G: false });
  const [tabValue, setTabValue] = useState("P");

  const [tempMethod, setTempMethod] = useState("FDR");
  const [tempThreshold, setTempThreshold] = useState("0.01");
  const [method, setMethod] = useState("FDR");
  const [threshold, setThreshold] = useState("0.01");

  // Evidence Modal 狀態
  const [detailModalOpen, setDetailModalOpen] = useState(false);
  const [selectedDetailRow, setSelectedDetailRow] = useState(null);
  const [evidenceRows, setEvidenceRows] = useState([]);
  const [evidenceLoading, setEvidenceLoading] = useState(false);

  const handleOpenDetail = async (row) => {
    setSelectedDetailRow(row);
    setDetailModalOpen(true);
    setEvidenceLoading(true);
    setEvidenceRows([]);

    try {
      const genes = inputText;
      const term_id = row["Term ID"];
      const data = await fetchTermEvidence({ term_id, genes });
      
      const formattedRows = (data.rows || []).map((r, idx) => ({
        id: `${r.WormBase_ID}-${r.Evidence_Code}-${idx}`,
        ...r
      }));
      setEvidenceRows(formattedRows);
    } catch (e) {
      console.error("Failed to fetch details:", e);
    } finally {
      setEvidenceLoading(false);
    }
  };

  const handleCloseDetail = () => {
    setDetailModalOpen(false);
    setSelectedDetailRow(null);
  };

  const handleTabChange = (key) => setTabValue(key);

  const handleAspectChange = (event) => {
    const { name, checked } = event.target;
    setSelectedAspects((prev) => ({ ...prev, [name]: checked }));
  };

  const handleAnalyze = async () => {
    if (!inputText.trim()) return;
    if (!Object.values(selectedAspects).some((v) => v)) {
      setError("Please select at least one feature.");
      return;
    }
    setLoading(true);
    setError(null);
    setRawResult(null);
    try {
      const runFpc = selectedAspects.F || selectedAspects.P || selectedAspects.C;
      const runGeneGroup = selectedAspects.G;
      const data = await analyzeWormBase(inputText, runFpc, runGeneGroup);
      setRawResult(data);

      const allOptions = [...ASPECTS, FGG_FEATURE];
      const firstAvailable = allOptions.find((a) => selectedAspects[a.key])?.key;
      if (firstAvailable) setTabValue(firstAvailable);
    } catch (err) {
      console.error(err);
      setError(err.message || "Analysis failed.");
    } finally {
      setLoading(false);
    }
  };

  const handleRefreshFilters = () => {
    setMethod(tempMethod);
    setThreshold(tempThreshold);
  };

  const handleSample = () => {
    setListName("Sample_Worm_Gene_List");
    setInputText(["WBGene00002324", "WBGene00003164", "WBGene00004343"].join("\n"));
  };

  const displayResult = useMemo(() => {
    if (!rawResult) return null;
    const threshVal = threshold === "None" ? 1.0 : parseFloat(threshold);

    const normalizeRow = (row, direction, isGeneGroup) => {
      let p_raw, p_fdr, p_bon, termId, termName, raw_fc, log2_fc, neg_log10_p;

      if (isGeneGroup) {
        p_raw = row["P-value"];
        p_fdr = row["FDR P-value"];
        p_bon = row["Bonferroni P-value"];
        termId = row["Feature Name"];
        termName = row["Feature Name"];
        raw_fc = row["Fold Change"];
        log2_fc = row["log2(Fold Change)"];
        neg_log10_p = row["-log10(P-value)"];
      } else {
        p_raw = row["p_raw"];
        p_fdr = row["p_fdr"];
        p_bon = row["p_bon"];
        termId = row["GO ID"];
        termName = row["GO ID"];
        log2_fc = row["Fold Change"] || 0;
        raw_fc = log2_fc === 0 ? 0 : Math.pow(2, log2_fc);
        neg_log10_p = p_raw > 0 ? -Math.log10(p_raw) : 300;
      }

      let display_p = p_raw;
      if (method === "FDR") display_p = p_fdr;
      if (method === "Bonferroni") display_p = p_bon;

      return {
        ...row,
        direction,
        "Term ID": termId,
        "Term Name": termName,
        p_raw, p_fdr, p_bon, display_p,
        "Fold Change": raw_fc, log2_fc, neg_log10_p,
        "Expected Ratio Text": isGeneGroup ? row["Expected Ratio"] : null,
        "Observed Ratio Text": isGeneGroup ? row["Observed Ratio"] : null,
      };
    };

    const processList = (list, direction, isGeneGroup) => {
      if (!list) return [];
      return list.map((row) => normalizeRow(row, direction, isGeneGroup)).filter((row) => row.display_p <= threshVal);
    };

    const mergeData = (aspectData, isGeneGroup) => {
      if (!aspectData) return [];
      const enriched = processList(aspectData.enriched, "Enriched", isGeneGroup);
      const depleted = processList(aspectData.depleted, "Depleted", isGeneGroup);
      return [...enriched, ...depleted].sort((a, b) => a.display_p - b.display_p);
    };

    const fpc = rawResult.fpc_results || {};
    const fgg = rawResult.gene_group_results || {};

    return {
      meta: rawResult.meta,
      F: mergeData(fpc.F, false),
      P: mergeData(fpc.P, false),
      C: mergeData(fpc.C, false),
      G: mergeData(fgg, true),
    };
  }, [rawResult, method, threshold]);

  const getCurrentData = () => {
    if (!displayResult) return [];
    if (tabValue === "F" && selectedAspects.F) return displayResult.F;
    if (tabValue === "P" && selectedAspects.P) return displayResult.P;
    if (tabValue === "C" && selectedAspects.C) return displayResult.C;
    if (tabValue === "G" && selectedAspects.G) return displayResult.G;
    return [];
  };

  const getFeatureLabel = (key) => {
    const found = [...ASPECTS, FGG_FEATURE].find((a) => a.key === key);
    return found ? found.btnLabel : "";
  };

  const handleDownload = () => {
    const data = getCurrentData();
    if (!data || data.length === 0) {
      alert("No data available to download.");
      return;
    }
    const csvRows = [
      ["Feature Name", "Observed Ratio", "Expected Ratio", "Fold Change", "P-value", "FDR P-value", "Bonferroni P-value", "log2(Fold Change)", "-log10(P-value)", "-log10(FDR P-value)"],
    ];
    data.forEach((row) => {
      const obsRatioStr = row["Observed Ratio Text"] || `${row["Input Hit"]}/${row["TotalInput"]} (${((row["Input Hit"] / row["TotalInput"]) * 100).toFixed(2)}%)`;
      const expRatioStr = row["Expected Ratio Text"] || `${row["Gene Count"]}/${row["TotalInAspect"]} (${((row["Gene Count"] / row["TotalInAspect"]) * 100).toFixed(2)}%)`;
      const featureName = row["Term ID"]?.startsWith("GO:") ? `${row["Term ID"]} (${row["Term Name"]})` : row["Term Name"];
      csvRows.push([`"${featureName}"`, `"${obsRatioStr}"`, `"${expRatioStr}"`, row["Fold Change"], row["p_raw"], row["p_fdr"], row["p_bon"], row["log2_fc"], row["neg_log10_p"], row["neg_log10_fdr"] || row["neg_log10_p"]]);
    });
    const csvContent = csvRows.map((e) => e.join(",")).join("\n");
    const blob = new Blob(["\uFEFF" + csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    const tabNameMap = { F: "MF", P: "BP", C: "CC", G: "FGG" };
    link.setAttribute("download", `${tabNameMap[tabValue] || "Export"}_Enrichment.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const handleDownloadEvidence = () => {
    if (evidenceRows.length === 0) return;
    const csvRows = [["WormBase ID", "Symbol", "GO ID", "Reference", "Evidence Code"]];
    evidenceRows.forEach((r) => {
      csvRows.push([r.WormBase_ID, r.Symbol, r.GO_ID, r.Reference, r.Evidence_Code]);
    });
    const csvContent = csvRows.map((e) => e.join(",")).join("\n");
    const blob = new Blob(["\uFEFF" + csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.setAttribute("download", `${selectedDetailRow["Term ID"]}_Evidence.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const styles = {
    sectionBorder: "1px solid #ddd", headerBg: "#f9f9f9", accentBlueBg: "#e3f2fd", accentBlueText: "#1976d2", btnGreen: "#1b5e20", outerBorder: "1px solid #ccc",
    mainHeaderBg: "#dcdcdc", subHeaderBg: "#f5f5f5", tableBorderColor: "#e0e0e0", labelColumnBg: "#f0f4f8", darkBarBg: "#1c384e",
    tabBtnStyle: (isActive) => ({
      textTransform: "none", boxShadow: "none", fontSize: "0.85rem", px: 2, border: "1px solid", borderRadius: 1,
      borderColor: isActive ? "#0d6efd" : "#bbb", bgcolor: isActive ? "#0d6efd" : "#e0e0e0", color: isActive ? "#fff" : "#333",
      "&:hover": { bgcolor: isActive ? "#0b5ed7" : "#d5d5d5", borderColor: isActive ? "#0b5ed7" : "#999" },
    }),
  };

  return (
    <Box sx={{ p: 4, maxWidth: 1400, margin: "0 auto" }}>
      <Grid container spacing={3}>
        {/* Step 1 */}
        <Grid item xs={12} md={5}>
          <Paper elevation={0} sx={{ border: styles.sectionBorder, height: "100%" }}>
            <Box sx={{ p: 2, borderBottom: styles.sectionBorder, bgcolor: styles.headerBg }}>
              <Typography variant="subtitle1" fontWeight="bold">Step1. Input gene list</Typography>
            </Box>
            <Box sx={{ p: 2 }}>
              <Box sx={{ display: "flex", gap: 1, mb: 2 }}>
                <TextField size="small" fullWidth placeholder="List Name" value={listName} onChange={(e) => setListName(e.target.value)} />
                <Button variant="contained" size="small" onClick={() => { setInputText(""); setListName(""); }} sx={{ bgcolor: "#757575" }}>Clear all</Button>
              </Box>
              <TextField multiline rows={10} fullWidth variant="outlined" value={inputText} onChange={(e) => setInputText(e.target.value)} sx={{ mb: 2, bgcolor: "#fff" }} placeholder="Paste IDs here..." />
              <Button variant="contained" color="primary" sx={{ textTransform: "none" }} onClick={handleSample}>Load Sample Input</Button>
            </Box>
          </Paper>
        </Grid>

        {/* Step 2 */}
        <Grid item xs={12} md={7}>
          <Paper elevation={0} sx={{ border: styles.sectionBorder, height: "100%" }}>
            <Box sx={{ p: 2, borderBottom: styles.sectionBorder, bgcolor: styles.headerBg }}>
              <Typography variant="subtitle1" fontWeight="bold">Step2. Select biological feature(s)</Typography>
            </Box>
            <Box sx={{ p: 2 }}>
              <Accordion defaultExpanded elevation={0} sx={{ border: "1px solid #bbdefb", mb: 2 }}>
                <AccordionSummary expandIcon={<ExpandMoreIcon sx={{ color: styles.accentBlueText }} />} sx={{ bgcolor: styles.accentBlueBg, minHeight: 48 }}>
                  <Typography color={styles.accentBlueText} fontWeight="bold">Ontology and Annotation Features</Typography>
                </AccordionSummary>
                <AccordionDetails sx={{ pt: 2 }}>
                  <Typography variant="body2" sx={{ textDecoration: "underline", color: "primary.main", mb: 1 }}>Gene Ontology (GO)</Typography>
                  <FormGroup row sx={{ ml: 1 }}>
                    {ASPECTS.map((aspect) => (
                      <FormControlLabel key={aspect.key} control={<Checkbox size="small" checked={selectedAspects[aspect.key]} onChange={handleAspectChange} name={aspect.key} />} label={<Typography variant="body2">{aspect.label}</Typography>} sx={{ width: "45%" }} />
                    ))}
                  </FormGroup>
                </AccordionDetails>
              </Accordion>

              <Accordion defaultExpanded elevation={0} sx={{ border: "1px solid #bbdefb", mb: 2 }}>
                <AccordionSummary expandIcon={<ExpandMoreIcon sx={{ color: styles.accentBlueText }} />} sx={{ bgcolor: styles.accentBlueBg, minHeight: 48 }}>
                  <Typography color={styles.accentBlueText} fontWeight="bold">Functional Gene Group Features</Typography>
                </AccordionSummary>
                <AccordionDetails sx={{ pt: 2 }}>
                  <FormGroup row sx={{ ml: 1 }}>
                    <FormControlLabel control={<Checkbox size="small" checked={selectedAspects.G} onChange={handleAspectChange} name="G" />} label={<Typography variant="body2">{FGG_FEATURE.label}</Typography>} />
                  </FormGroup>
                </AccordionDetails>
              </Accordion>
            </Box>
          </Paper>
        </Grid>
      </Grid>

      {/* 執行按鈕 */}
      <Box sx={{ mt: 4, display: "flex", justifyContent: "center" }}>
        <Button variant="contained" size="large" onClick={handleAnalyze} disabled={loading} startIcon={!loading && <SendIcon />} sx={{ minWidth: 250, bgcolor: styles.btnGreen, fontSize: "1.1rem", textTransform: "none", py: 1.5, "&:hover": { bgcolor: "#144a17" } }}>
          {loading ? <CircularProgress size={26} color="inherit" /> : "Run Enrichment Analysis"}
        </Button>
      </Box>

      {error && <Alert severity="error" sx={{ mt: 2, mx: "auto", maxWidth: 800 }}>{error}</Alert>}

      {/* ======================= 分析結果區塊 ======================= */}
      {displayResult && (
        <Box sx={{ mt: 5 }}>
          <Paper elevation={0} sx={{ border: styles.outerBorder, mb: 4, overflow: "hidden" }}>
            <Box sx={{ bgcolor: styles.subHeaderBg, p: 1.5, display: "flex", alignItems: "center", borderBottom: styles.tableBorderColor }}>
              <ManageAccountsIcon sx={{ mr: 1, fontSize: 28 }} />
              <Typography variant="h6" fontWeight="bold">User's Specification</Typography>
            </Box>

            <Box>
              <Box sx={{ display: "flex", borderBottom: `1px solid ${styles.tableBorderColor}` }}>
                <Box sx={{ flex: "0 0 35%", bgcolor: styles.labelColumnBg, p: 2, textAlign: "center", borderRight: `1px solid ${styles.tableBorderColor}` }}>
                  <Typography variant="body2">Name of each input list</Typography>
                </Box>
                <Box sx={{ flex: 1, p: 2, textAlign: "center" }}><Typography variant="body2">{listName || "N/A"}</Typography></Box>
              </Box>
              <Box sx={{ display: "flex", borderBottom: `1px solid ${styles.tableBorderColor}` }}>
                <Box sx={{ flex: "0 0 35%", bgcolor: styles.labelColumnBg, p: 2, textAlign: "center", borderRight: `1px solid ${styles.tableBorderColor}` }}>
                  <Typography variant="body2"># of genes in each input list</Typography>
                </Box>
                <Box sx={{ flex: 1, p: 2, textAlign: "center" }}><Typography variant="body2">{displayResult.meta.mapped_live}</Typography></Box>
              </Box>

              <Box sx={{ bgcolor: styles.darkBarBg, p: 1, textAlign: "center" }}><Typography variant="body2" color="white">See the analysis result of a chosen feature</Typography></Box>

              <Box sx={{ display: "flex", borderBottom: `1px solid ${styles.tableBorderColor}` }}>
                <Box sx={{ flex: "0 0 35%", p: 2, display: "flex", alignItems: "center", justifyContent: "center", borderRight: `1px solid ${styles.tableBorderColor}` }}>
                  <Typography variant="body2">Ontology and Annotation Features</Typography>
                </Box>
                <Box sx={{ flex: 1, p: 2, display: "flex", gap: 2, flexWrap: "wrap" }}>
                  {ASPECTS.map((aspect) => selectedAspects[aspect.key] && (
                    <Button key={aspect.key} onClick={() => handleTabChange(aspect.key)} sx={styles.tabBtnStyle(tabValue === aspect.key)}>{aspect.btnLabel}</Button>
                  ))}
                </Box>
              </Box>

              <Box sx={{ display: "flex", borderBottom: `1px solid ${styles.tableBorderColor}` }}>
                <Box sx={{ flex: "0 0 35%", p: 2, display: "flex", alignItems: "center", justifyContent: "center", borderRight: `1px solid ${styles.tableBorderColor}` }}>
                  <Typography variant="body2">Functional Gene Group Features</Typography>
                </Box>
                <Box sx={{ flex: 1, p: 2, display: "flex", gap: 2, flexWrap: "wrap" }}>
                  {selectedAspects.G && <Button onClick={() => handleTabChange("G")} sx={styles.tabBtnStyle(tabValue === "G")}>{FGG_FEATURE.btnLabel}</Button>}
                </Box>
              </Box>
            </Box>
          </Paper>

          <Box sx={{ mb: 2, display: "flex", alignItems: "center" }}>
            <CheckIcon sx={{ mr: 1, fontSize: 32 }} />
            <Typography variant="h5">Results of the analyzed feature: <b>{getFeatureLabel(tabValue)}</b></Typography>
          </Box>

          <Paper elevation={0} sx={{ border: styles.outerBorder, mb: 3, bgcolor: "#fbfbfb" }}>
            <Box sx={{ p: 1.5, borderBottom: styles.tableBorderColor, display: "flex", alignItems: "center" }}>
              <SettingsIcon sx={{ mr: 1, fontSize: 20 }} />
              <Typography variant="body2" fontWeight="bold">Specify the correction method and p-value cutoff for multiple hypotheses testing</Typography>
            </Box>
            <Box sx={{ p: 3, display: "flex", gap: 5, alignItems: "flex-end" }}>
              <Box>
                <Typography variant="body2" fontWeight="bold" mb={1}>Correction method:</Typography>
                <Select size="small" value={tempMethod} onChange={(e) => setTempMethod(e.target.value)} sx={{ bgcolor: "white", minWidth: 200 }}>
                  <MenuItem value="None">None</MenuItem><MenuItem value="FDR">FDR</MenuItem><MenuItem value="Bonferroni">Bonferroni</MenuItem>
                </Select>
              </Box>
              <Box>
                <Typography variant="body2" fontWeight="bold" mb={1}>Corrected p-value cutoff:</Typography>
                <Select size="small" value={tempThreshold} onChange={(e) => setTempThreshold(e.target.value)} sx={{ bgcolor: "white", minWidth: 200 }}>
                  <MenuItem value="0.05">0.05</MenuItem><MenuItem value="0.01">0.01</MenuItem><MenuItem value="None">No Cutoff</MenuItem>
                </Select>
              </Box>
              <Box>
                <Button variant="contained" onClick={handleRefreshFilters} startIcon={<RefreshIcon />} sx={{ bgcolor: "#00bcd4", color: "black", textTransform: "none", fontWeight: "bold", "&:hover": { bgcolor: "#00acc1" } }}>Refresh</Button>
              </Box>
            </Box>
          </Paper>

          <Paper elevation={0} sx={{ border: styles.outerBorder, p: 2 }}>
            <Box sx={{ mb: 2, p: 2, bgcolor: "#fafafa", borderRadius: 1 }}>
              <Typography variant="body2"><strong>A:</strong> # of input genes with this term</Typography>
              <Typography variant="body2"><strong>B:</strong> # of input genes</Typography>
              <Typography variant="body2"><strong>C:</strong> # of genes (in the <em>C. elegans</em> genome) with this term</Typography>
              <Typography variant="body2"><strong>D:</strong> # of genes in the <em>C. elegans</em> genome ({displayResult.meta.total_background} genes)</Typography>
              <FoldEnrichmentFormula />
            </Box>
            <Box sx={{ mb: 2 }}>
              <Button variant="outlined" startIcon={<DownloadIcon />} size="small" onClick={handleDownload} sx={{ textTransform: "none", color: "#333", borderColor: "#ccc" }}>Download as a .csv file</Button>
            </Box>
            <AspectResultPanel data={getCurrentData()} methodLabel={method === "None" ? "Raw" : method} onObservedClick={handleOpenDetail} />
          </Paper>

          <Paper elevation={0} sx={{ border: styles.outerBorder, mt: 3 }}>
            <Box sx={{ bgcolor: styles.subHeaderBg, p: 1, display: "flex", alignItems: "center", borderBottom: styles.tableBorderColor }}>
              <BarChartIcon sx={{ mr: 1 }} />
              <Typography variant="subtitle1" fontWeight="bold">Graphic View</Typography>
            </Box>
            <GraphicViewPanel data={getCurrentData()} threshold={threshold} />
          </Paper>
        </Box>
      )}

      {/* ======================= Evidence 詳細資訊 Modal ======================= */}
      <Dialog
        open={detailModalOpen}
        onClose={handleCloseDetail}
        maxWidth="lg"
        fullWidth
        scroll="paper"
        PaperProps={{ sx: { bgcolor: "#f5f5f5", maxHeight: "90vh" } }}
      >
        <DialogTitle sx={{ display: "flex", justifyContent: "space-between", alignItems: "center", bgcolor: "#fff", borderBottom: "1px solid #ddd", py: 1.5 }}>
          <Typography variant="h6" component="div" fontWeight="bold">Enrichment Analysis Result</Typography>
          <IconButton onClick={handleCloseDetail} size="small" sx={{ border: "1px solid #ccc", borderRadius: 1 }}>
            <CloseIcon />
          </IconButton>
        </DialogTitle>

        <DialogContent sx={{ p: 3, bgcolor: "#f9f9f9" }}>
          {selectedDetailRow && displayResult && (
            <Box sx={{ display: "flex", flexDirection: "column", gap: 3 }}>
              
              {/* 區塊 1: User's Specification */}
              <Paper variant="outlined" sx={{ bgcolor: "#fff", borderRadius: 2 }}>
                <Box sx={{ bgcolor: "#f5f5f5", p: 1.5, borderBottom: "1px solid #e0e0e0", display: "flex", alignItems: "center" }}>
                  <ManageAccountsIcon sx={{ mr: 1, fontSize: 20, color: "#5e35b1" }} />
                  <Typography variant="subtitle2" fontWeight="bold">User's Specification</Typography>
                </Box>
                <Table size="small">
                  <TableBody>
                    <TableRow>
                      <TableCell sx={{ width: "30%", fontWeight: "bold", borderRight: "1px solid #e0e0e0" }}>Analyzed Term</TableCell>
                      <TableCell sx={{ color: "#1976d2", textDecoration: "underline", cursor: "pointer" }}>{selectedDetailRow["Term ID"]}</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ fontWeight: "bold", borderRight: "1px solid #e0e0e0" }}>Analyzed List</TableCell>
                      <TableCell>{listName || "User Input List"}</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ fontWeight: "bold", borderRight: "1px solid #e0e0e0" }}># of genes in the list</TableCell>
                      <TableCell>{displayResult.meta.mapped_live}</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell sx={{ fontWeight: "bold", borderRight: "1px solid #e0e0e0" }}>Multiple hypothesis testing</TableCell>
                      <TableCell>{method.toUpperCase()} correction, Cutoff={threshold}</TableCell>
                    </TableRow>
                  </TableBody>
                </Table>
              </Paper>

              {/* 區塊 2: Analysis Result Summary */}
              <Paper variant="outlined" sx={{ bgcolor: "#fff", borderRadius: 2 }}>
                <Box sx={{ bgcolor: "#f5f5f5", p: 1.5, borderBottom: "1px solid #e0e0e0", display: "flex", alignItems: "center" }}>
                  <PushPinIcon sx={{ mr: 1, fontSize: 20, color: "#d81b60" }} />
                  <Typography variant="subtitle2" fontWeight="bold">Analysis Result</Typography>
                </Box>
                <Table size="small">
                  <TableHead>
                    <TableRow>
                      <TableCell><b>Enriched Term</b></TableCell>
                      <TableCell><b>Observed Ratio</b></TableCell>
                      <TableCell><b>Expected Ratio</b></TableCell>
                      <TableCell><b>Trend</b></TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    <TableRow>
                      <TableCell sx={{ color: "#1976d2", textDecoration: "underline" }}>{selectedDetailRow["Term ID"]}</TableCell>
                      <TableCell>{selectedDetailRow["Input Hit"]}/{selectedDetailRow["TotalInput"]}</TableCell>
                      <TableCell>{selectedDetailRow["Expected Ratio Text"] || `${selectedDetailRow["Gene Count"]}/${selectedDetailRow["TotalInAspect"]} (${((selectedDetailRow["Gene Count"] / selectedDetailRow["TotalInAspect"]) * 100).toFixed(3)}%)`}</TableCell>
                      <TableCell>{selectedDetailRow.direction}</TableCell>
                    </TableRow>
                  </TableBody>
                </Table>
              </Paper>

              {/* 區塊 3: Evidence 的紅字標題 與 表格 */}
              <Box sx={{ mt: 1 }}>
                <Typography variant="h6" sx={{ color: "#c62828", fontWeight: "bold", mb: 2, textAlign: "left" }}>
                  {evidenceRows.length} evidence of the {selectedDetailRow["Input Hit"]} input genes that are associated with the GO term <Box component="span" sx={{ color: "#1976d2", textDecoration: "underline" }}>[{selectedDetailRow["Term ID"]}]</Box>
                </Typography>
                
                <Box sx={{ display: "flex", justifyContent: "center", mb: 2 }}>
                  <Button 
                    variant="contained" 
                    startIcon={<DownloadIcon />} 
                    onClick={handleDownloadEvidence}
                    sx={{ bgcolor: "#455a64", color: "white", textTransform: "none", fontWeight: "bold", borderRadius: 2, '&:hover': { bgcolor: "#37474f" } }}
                  >
                    Download (Evidence)
                  </Button>
                </Box>

                <Paper variant="outlined" sx={{ bgcolor: "#fff", height: 400, width: "100%" }}>
                  <DataGrid
                    rows={evidenceRows}
                    columns={evidenceColumns}
                    loading={evidenceLoading}
                    pageSizeOptions={[10, 25, 50]}
                    initialState={{ pagination: { paginationModel: { pageSize: 10 } } }}
                    disableRowSelectionOnClick
                    density="compact"
                    sx={{ border: "none" }}
                  />
                </Paper>
              </Box>

            </Box>
          )}
        </DialogContent>
      </Dialog>
    </Box>
  );
}