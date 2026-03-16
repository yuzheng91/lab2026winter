import React, { useState } from "react";
import {
  Box,
  Button,
  TextField,
  Typography,
  Paper,
  CircularProgress,
  Alert,
  Grid,
  FormGroup,
  FormControlLabel,
  Checkbox,
  Radio,
  RadioGroup,
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
  ToggleButton,
  ToggleButtonGroup,
} from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import DescriptionIcon from "@mui/icons-material/Description";
import TableViewIcon from "@mui/icons-material/TableView";
import InsertChartIcon from "@mui/icons-material/InsertChart";
import Plot from "react-plotly.js";

import { runWqfcPipeline, fetchWqfcGeneDetail } from "../api/wqfcApi";

const SAMPLE_1 = `YKL144C
YML060W
YKL110C
YIL096C
YNL227C
YLR336C
YGR251W
YGR054W
YHR042W
YPR034W
YLL036C
YML082W
YOR233W
YMR277W
YOR217W
YLL038C
YLR433C
YBL035C
YOR188W
YNL102W
YDR097C
YNL302C
YOR063W
YDR211W
YKL154W
YDR190C
YGR229C
YIL021W
YNL153C
YHR144C
YAL035W
YKR079C
YDR449C
YDR161W
YJR003C
YML018C
YLR405W
YGR195W
YGL070C
YJR043C
YIL110W
YNR009W
YBR089W
YLR183C
YER070W
YBR030W
YNL228W
YMR269W
YNR054C
YLR400W
YDL151C
YNL303W
YKR025W
YDR152W
YNL292W
YOL097C
YLR340W
YBR106W
YOR253W
YBL077W
YHR013C
YDR233C
YPR080W
YBR118W
YLR134W
YLR044C
YFL045C
YER009W
YHR133C
YLR447C
YHR026W
YKL080W
YOR209C
YMR079W
YAL038W
YBR249C
YPR060C
YGR151C
YGR152C
YGR228W
YGR081C
YIL079C
YIL104C
YDR440W
YNR024W
YOR252W
YKR060W
YHL013C
YOL124C
YPL043W
YAL025C
YBR266C
YLR017W
YOL039W
YOL120C
YPL207W
YOR096W
YBR252W
YGR078C
YGL148W
YOR271C
YKR043C
YOL040C
YOL022C
YNR003C
YLR073C
YLR146C
YPR016C
YNL255C
YOR146W
YIL019W
YHR148W
YLR074C
YPL044C
YPR137W
YDL060W
YAL059W
YER118C
YGL097W
YNL066W
YHR025W
YML106W
YGL008C
YEL042W
YBR283C
YCR034W
YMR215W
YDL084W
YBR154C
YDR023W
YDR429C
YLR150W
YNL010W
YPL037C
YHR193C
YJL208C
YOR277C
YGR285C
YHL011C
YPR041W
YPR136C
YPL217C
YMR128W
YPR112C
YLL034C
YPL263C
YKL205W
YPR010C
YJL050W
YJL080C
YOR361C
YGR265W
YGR162W
YBR079C
YOR287C
YER082C
YFL023W
YLR051C
YMR310C
YGR283C
YMR014W
YBR155W
YDL150W
YER127W
YMR239C
YGR280C
YHR088W
YNL186W
YBR267W
YNL299W
YDL062W
YKL143W
YDR120C
YGR173W
YML043C
YLR243W
YKR026C
YPL245W
YDL111C
YDR300C
YER049W
YNL163C
YJL198W
YER036C
YEL040W
YER043C
YAR075W
YAR073W
YLR432W
YHR216W
YML056C
YMR309C
YPR187W
YJR063W
YMR235C
YPR118W
YJR007W
YLR372W
YOR260W
YLR291C
YLR293C
YLR083C
YCR053W
YHR064C
YDR037W
YPL160W
YGR094W
YNL247W
YBR121C
YIL078W
YHR019C
YLR060W
YPR033C
YMR307W
YPR074C
YCL045C
YJR143C
YGR185C
YLR167W
YHR128W
YJL191W
YMR146C
YDL081C
YJL136C
YML022W
YLR367W
YOR167C
YMR121C
YLR029C
YMR116C
YDR454C
YJL138C
YLR264W
YKL180W
YKR059W
YOR234C
YDL083C
YLR325C
YMR194W
YMR230W
YOR293W
YBL072C
YPL143W
YML026C
YGR085C
YBL092W
YBL027W
YBR084C-A
YBR181C
YGL189C
YGR027C
YER131W
YGR118W
YKR094C
YPR132W
YIL052C
YPL090C
YER056C-A
YNL162W
YGR034W
YML063W
YNL209W
YLL045C
YGL030W
YML073C
YPL142C
YKL006W
YNL119W
YLR448W
YLR388W
YHL033C
YER074W
YBR189W
YBR191W
YGL135W
YPL220W
YLR048W
YLL044W
YGL123W
YIL018W
YFR031C-A
YDR025W
YJR123W
YDR417C
YLR062C
YNL301C
YPL198W
YGL076C
YJL177W
YER102W
YPR102C
YMR242C
YGL147C
YHR203C
YBR048W
YNL178W
YGR148C
YGL031C
YGR214W
YIL133C
YOR369C
YDL082W
YEL054C
YDR064W
YLR344W
YPR044C
YGL102C
YNL069C
YJL190C
YML024W
YMR142C
YPL079W
YPL081W
YNL096C
YJR145C
YLR061W
YHR010W
YDR450W
YHR141C
YPR043W
YER117W
YBL087C
YHL001W
YKR057W
YDR471W
YKL056C
YLR441C
YOR312C
YDR447C
YJR094W-A
YLR185W
YDR382W
YDR418W
YDR012W
YCR031C
YBR031W
YEL026W
YLR075W
YIL148W
YLR339C
YKL081W
YDL191W
YDL061C
YDR500C
YDL075W
YER025W
YDR321W
YLR397C
YKR092C
YJR071W
YNL075W
YNL114C
YNL113W
YNL022C
YLR186W
YDR184C
YPL211W
YFL002C
YPR144C
YHR196W
YIR012W
YLR401C
YMR093W
YDL167C
YDR365C
YCL054W
YNL061W
YOR294W
YOR145C
YHR169W
YJL122W
YNL141W
YER110C
YJR041C
YDR341C
YOR168W
YMR308C
YDR385W
YGR264C
YPL226W
YDR091C
YLR249W
YGL120C
YMR217W
YPL183C
YNL313C
YJL109C
YCL037C
YMR229C
YBR084W
YHR020W
YMR049C
YDR412W
YDR312W
YKL078W
YOR243C
YDR465C
YFR001W
YDR361C
YKL099C
YBL028C
YIL127C
YHR066W
YHR197W
YOR272W
YJR070C
YLR222C
YGR123C
YKL181W
YCL059C
YKR081C
YKL021C
YNL110C
YLR221C
YDR413C
YJL148W
YCR016W
YLR003C
YHR052W
YKL009W
YNL002C
YDR101C
YDR324C
YDR083W
YGL099W
YGL111W
YDL051W
YMR131C
YPR110C
YMR290C
YBL039C
YML093W
YGL078C
YLL008W
YER006W
YDR496C
YGR145W
YGR245C
YNL175C
YNL174W
YLR196W
YOR206W
YGR103W
YLR276C
YGR160W
YNL182C
YLR409C
YLR129W
YLR002C
YNL132W
YNL248C
YGR128C
YLR009W
YNL062C
YGL171W
YKR056W
YBR247C
YJR002W
YLR198C
YOR309C
YLR197W
YLR175W
YOR310C
YLR449W
YHR170W
YBL024W
YBR143C
YDR165W
YDR398W
YPR190C
YCR072C
YJL069C
YDL201W
YHR065C
YLL011W
YER126C
YDR087C
YDL153C
YPR142C
YPR143W
YPL146C
YKL172W
YNL308C
YKL082C
YGR271C-A
YGL029W
YPL266W
YOL077C
YBR142W
YBR034C
YLR435W
YGR187C
YHR089C
YCR057C
YBL068W
YPR163C
YOR224C
YLR413W
YJR124C
YPL212C
YPL086C
YNL256W
YIR026C
YGR200C
YJL125C
YNL120C
YKR024C
YBL054W
YNL124W
YBR061C
YLR063W
YML080W
YGL169W
YDL063C
YGR158C
YIL020C
YHR143W-A
YGR083C
YER165W
YAL036C
YDL031W
YJL188C
YJL189W
YDL136W
YNL067W
YPL238C
YPL237W
YOR276W
YNR043W
YPR145W
YGL213C
YPL273W
YMR321C
YMR143W
YJL183W
YHR068W
YJR064W
YDL143W
YJL014W
YJL002C
YGR124W
YMR212C
YHL039W
YMR243C
YBR238C
YDR346C
YDR280W
YCL031C
YBR025C
YOL041C
YLL035W
YOR004W
YOR001W
YOR210W
YAL003W
YJL033W
YGL103W
YLR076C
YNL151C
YLR172C
YFL022C
YDR172W
YOR169C
YOL121C
YLR333C
YOL127W
YOR021C
YOR182C
YHR021C
YKL156W
YDL050C
YOR095C
YLR065C
YBR246W
YGL225W
YNR053C
YPL235W
YNL001W
YDR339C
YJR014W
YLR412W
YKL122C
YMR234W
YCR035C
YPR169W`;

const SAMPLE_2 = `YKL121W
YLR001C
YNL183C
YNL241C
YNL208W
YKR066C
YPR067W
YKR067W
YKL100C
YMR304W
YKL193C
YDR032C
YCR004C
YML004C
YKL065C
YLR290C
YKL067W
YGR209C
YDR453C
YML131W
YNL134C
YOL151W
YOL150C
YKR011C
YKL103C
YDL124W
YMR173W-A
YMR173W
YDR513W
YMR315W
YKL142W
YJR025C
YGR127W
YLR250W
YDR368W
YML110C
YLR370C
YJR104C
YDR512C
YGL010W
YEL005C
YDR286C
YIL120W
YMR155W
YDR350C
YMR253C
YMR139W
YMR169C
YOR052C
YOL048C
YER142C
YHR112C
YGL087C
YMR197C
YDL021W
YLL039C
YOR185C
YOR173W
YGL248W
YPR026W
YOR317W
YMR136W
YJL141C
YKR058W
YOL082W
YNL305C
YBR287W
YER141W
YBR183W
YDR255C
YIL097W
YHR171W
YKL124W
YMR041C
YDR254W
YGL156W
YBR285W
YGR008C
YKL151C
YDL110C
YEL060C
YHR087W
YGR043C
YDR272W
YDR085C
YNR001C
YBR139W
YIL124W
YBR026C
YNL055C
YLR299W
YOR152C
YMR181C
YLR327C
YOR374W
YJL164C
YIL107C
YMR170C
YEL011W
YLR345W
YPL004C
YBR149W
YCL035C
YNL015W
YLR251W
YHR104W
YLL026W
YOL083W
YHR138C
YIR037W
YOR285W
YPL087W
YPR098C
YPL154C
YMR297W
YKL150W
YMR110C
YLR252W
YKL091C
YLR270W
YDR358W
YBR230C
YMR196W
YDR001C
YML128C
YFL014W
YNL274C
YLR178C
YJR096W
YGR248W
YPL230W
YLR149C
YGR088W
YDL204W
YBR053C
YBR056W
YBR126C
YCL042W
YCL040W
YGL037C
YML100W
YBR169C
YLR258W
YMR105C
YFR053C
YMR250W
YHL021C
YDR171W
YBL064C
YKR076W
YMR090W
YIR038C
YIR039C
YER079W
YER053C
YOR289W
YLR080W
YGR019W
YMR251W-A
YBR006W
YPL165C
YNL160W
YER037W
YOR220W
YHR097C
YJL057C
YBR269C
YMR081C
YKL093W
YHR096C
YBR072W
YIL136W
YJL048C
YDR059C
YPL054W
YLR142W
YIL101C
YJL066C
YNL115C
YMR174C
YDR070C
YMR271C
YGR201C
YJR008W
YPL186C
YJL103C
YLR312C
YKL026C
YPR150W
YPR151C
YML042W
YBL049W
YNL009W
YKR009C
YIL045W
YGL121C
YLL020C
YJL163C
YDL199C
YGL208W
YER150W
YNL200C
YKL133C
YJL161W
YDR533C
YPL166W
YPL203W
YLL019C
YNR002C
YDR096W
YCR091W
YER035W
YMR302C
YBL078C
YGL180W
YLR267W
YKL188C
YGR256W
YJL016W
YJL015C
YJL020C
YGR149W
YLR152C
YGR194C
YGR250C
YDR258C
YPL196W
YGL104C
YDR169C
YEL012W
YGL096W
YMR251W
YLR356W
YGL059W
YPL152W
YDL027C
YKR049C
YMR053C
YPR184W
YPR149W
YBR052C
YDR516C
YOR161C
YKL035W
YJR059W
YDR074W
YDL023C
YDL022W
YGR237C
YGL006W
YMR291W
YIL099W
YDR380W
YHR137W
YDL222C
YJR066W
YER103W
YBL075C
YNL194C
YAL061W
YIL113W
YNL195C
YOL084W
YGR023W
YGR161C
YJL144W
YJL142C
YER158C
YJL155C
YBR204C
YEL039C
YDR003W
YJL132W
YER038C
YDR475C
YDR287W
YFL044C
YDR204W
YCR061W
YER119C
YJL068C`;

const FEATURE_GROUPS = [
  {
    category: "Gene Features",
    features: [
      "CDS length",
      "5'UTR length",
      "3'UTR length",
      "# of pathways",
      "# of GO terms",
    ],
  },
  {
    category: "mRNA Features",
    features: ["Transcriptional plasticity", "mRNA half-life"],
  },
  {
    category: "Protein Features",
    features: ["# of PTMs", "Aromaticity score"],
  },
  { category: "Network Features", features: ["EPA network", "LEA network"] },
];

export default function WqfcAnalysis() {
  const [list1Text, setList1Text] = useState("");
  const [list2Text, setList2Text] = useState("");
  const [selectedFeatures, setSelectedFeatures] = useState({});

  // Step 3 狀態
  const [method, setMethod] = useState("Bonferroni");
  const [cutoffStr, setCutoffStr] = useState("10^-3");

  // 結果狀態
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [rawResult, setRawResult] = useState(null);

  // 視圖控制
  const [activeFeature, setActiveFeature] = useState("");
  const [viewMode, setViewMode] = useState("table"); // 'table' or 'figure'

  const handleFeatureChange = (event) => {
    setSelectedFeatures((prev) => ({
      ...prev,
      [event.target.name]: event.target.checked,
    }));
  };

  const handleLoadSample1 = () => {
    setList1Text(SAMPLE_1);
    // 預設勾選一個特徵，例如 5'UTR length，讓使用者一載入就能直接 Submit
    setSelectedFeatures((prev) => ({ ...prev, "5'UTR length": true }));
  };

  const handleLoadSample2 = () => {
    setList2Text(SAMPLE_2);
    // 預設勾選一個特徵，例如 5'UTR length，讓使用者一載入就能直接 Submit
    setSelectedFeatures((prev) => ({ ...prev, "5'UTR length": true }));
  };

  const handleAnalyze = async () => {
    const featuresToRun = Object.keys(selectedFeatures).filter(
      (k) => selectedFeatures[k],
    );
    if (!list1Text.trim() || !list2Text.trim() || featuresToRun.length === 0) {
      setError("Please input two gene lists and select at least one feature.");
      return;
    }

    setLoading(true);
    setError(null);
    setRawResult(null);
    try {
      let numericCutoff = 0.05;
      try {
        if (cutoffStr.includes("^")) {
          const parts = cutoffStr.split("^");
          numericCutoff = Math.pow(parseFloat(parts[0]), parseFloat(parts[1]));
        } else {
          numericCutoff = parseFloat(cutoffStr);
        }
      } catch (e) {}

      const data = await runWqfcPipeline(
        list1Text,
        list2Text,
        featuresToRun,
        method,
        numericCutoff,
      );
      if (data.status === "error") throw new Error(data.message);

      setRawResult(data);
      setActiveFeature(featuresToRun[0]);
    } catch (err) {
      setError(err.message || "Analysis failed.");
    } finally {
      setLoading(false);
    }
  };

  const getPValue = (testRow, testType) => {
    if (method === "Bonferroni") return testRow[`${testType}_bonf`];
    if (method === "FDR") return testRow[`${testType}_fdr`];
    return testRow[`${testType}_test`];
  };

  const formatPValue = (val) => (val != null ? val.toExponential(3) : "N/A");

  // CDF 計算函式 (補回)
  const generateCdfData = (arr) => {
    const sorted = [...arr].sort((a, b) => a - b);
    return {
      x: sorted,
      y: sorted.map((_, i) => (i + 1) / sorted.length),
    };
  };

  const renderResultSection = () => {
    if (!rawResult) return null;

    const cutoffValue = rawResult.summary_global.cutoff;
    const currentData = rawResult.results[activeFeature];
    const currentPlot = rawResult.plot_data[activeFeature];

    return (
      <Box sx={{ mt: 5 }}>
        <Typography
          variant="h5"
          color="#1976d2"
          sx={{ display: "flex", alignItems: "center", mb: 2 }}
        >
          <DescriptionIcon sx={{ mr: 1, fontSize: 30 }} /> Comparison Results
        </Typography>

        <TableContainer
          component={Paper}
          variant="outlined"
          sx={{ mb: 4, borderRadius: 1 }}
        >
          <Table size="small" sx={{ borderCollapse: "collapse" }}>
            <TableHead sx={{ bgcolor: "#f5f5f5" }}>
              <TableRow>
                <TableCell
                  colSpan={4}
                  sx={{ py: 1.5, borderBottom: "1px solid #ddd" }}
                >
                  <Typography variant="subtitle2" fontWeight="bold">
                    ➔ User's Specification
                  </Typography>
                </TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {/* 基本資訊列 */}
              <TableRow>
                <TableCell
                  width="20%"
                  sx={{
                    fontWeight: "bold",
                    borderRight: "1px solid #ddd",
                    borderBottom: "1px solid #ddd",
                  }}
                >
                  # of genes in L1
                </TableCell>
                <TableCell
                  width="30%"
                  sx={{
                    borderRight: "1px solid #ddd",
                    borderBottom: "1px solid #ddd",
                  }}
                >
                  {rawResult.summary_global.L1_total}
                </TableCell>
                <TableCell
                  width="20%"
                  sx={{
                    fontWeight: "bold",
                    borderRight: "1px solid #ddd",
                    borderBottom: "1px solid #ddd",
                  }}
                >
                  # of genes in L2
                </TableCell>
                <TableCell width="30%" sx={{ borderBottom: "1px solid #ddd" }}>
                  {rawResult.summary_global.L2_total}
                </TableCell>
              </TableRow>
              <TableRow>
                <TableCell
                  sx={{ fontWeight: "bold", borderRight: "1px solid #ddd" }}
                >
                  Multiple hypotheses testing
                </TableCell>
                <TableCell colSpan={3}>
                  {method === "None" ? "No correction" : `${method} correction`}
                  : p-value cutoff = {cutoffStr}
                </TableCell>
              </TableRow>

              {/* ⭐ 滿版深藍色標題列 (利用 colSpan={4} 解決擠在左邊的問題) */}
              <TableRow>
                <TableCell
                  colSpan={4}
                  sx={{
                    bgcolor: "#1a4659",
                    color: "white",
                    textAlign: "center",
                    fontWeight: "bold",
                    py: 1.5,
                    borderTop: "2px solid #fff",
                  }}
                >
                  See the testing result of a chosen quantitative feature
                </TableCell>
              </TableRow>

              {/* 特徵選擇列 */}
              {FEATURE_GROUPS.map((group) => {
                const groupSelectedFeatures = group.features.filter(
                  (f) => selectedFeatures[f],
                );
                if (groupSelectedFeatures.length === 0) return null;
                return (
                  <TableRow key={group.category}>
                    <TableCell
                      sx={{
                        borderRight: "1px solid #ddd",
                        borderBottom: "1px solid #ddd",
                        bgcolor: "#fdfdfd",
                      }}
                    >
                      {group.category}
                    </TableCell>
                    <TableCell
                      colSpan={3}
                      sx={{ borderBottom: "1px solid #ddd" }}
                    >
                      {groupSelectedFeatures.map((feat) => (
                        <Button
                          key={feat}
                          variant={
                            activeFeature === feat ? "contained" : "outlined"
                          }
                          size="small"
                          onClick={() => setActiveFeature(feat)}
                          sx={{
                            m: 0.5,
                            textTransform: "none",
                            borderRadius: 1,
                          }}
                        >
                          {feat}
                        </Button>
                      ))}
                    </TableCell>
                  </TableRow>
                );
              })}
            </TableBody>
          </Table>
        </TableContainer>

        {/* --- 測試結果區域 --- */}
        <Box
          sx={{
            border: "1px solid #ddd",
            p: 2,
            bgcolor: "#fff",
            borderRadius: 1,
          }}
        >
          <Box sx={{ borderBottom: "1px solid #eee", pb: 2, mb: 2 }}>
            <Typography
              variant="subtitle2"
              fontWeight="bold"
              sx={{ bgcolor: "#f1f1f1", p: 1.5, borderRadius: 1 }}
            >
              ✔ Analyzed quantitative feature (QF): {activeFeature}
            </Typography>
          </Box>

          <ToggleButtonGroup
            value={viewMode}
            exclusive
            onChange={(e, newVal) => {
              if (newVal) setViewMode(newVal);
            }}
            size="small"
            sx={{ mb: 3 }}
          >
            <ToggleButton value="table">
              <TableViewIcon sx={{ mr: 1 }} /> Table View
            </ToggleButton>
            <ToggleButton value="figure">
              <InsertChartIcon sx={{ mr: 1 }} /> Figure View
            </ToggleButton>
          </ToggleButtonGroup>

          {/* 表格視圖 */}
          {viewMode === "table" && currentData && (
            <Box>
              {/* 略過重複的 Table View 程式碼，維持你之前的樣式即可，這部分原本沒問題 */}
              <Typography
                variant="subtitle2"
                sx={{
                  bgcolor: "#f5f5f5",
                  p: 1,
                  borderTop: "1px solid #ccc",
                  borderLeft: "1px solid #ccc",
                  borderRight: "1px solid #ccc",
                }}
              >
                📊 Summary
              </Typography>
              <TableContainer
                component={Paper}
                variant="outlined"
                sx={{ borderRadius: 0, borderTop: 0, mb: 4 }}
              >
                <Table size="small">
                  <TableHead>
                    <TableRow>
                      <TableCell width="40%"></TableCell>
                      <TableCell width="30%">
                        <strong>First Gene List (L1)</strong>
                      </TableCell>
                      <TableCell width="30%">
                        <strong>Second Gene List (L2)</strong>
                      </TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    <TableRow>
                      <TableCell>
                        <strong># of genes which have this feature data</strong>
                      </TableCell>
                      <TableCell>
                        <a href="#" style={{ color: "#1976d2" }}>
                          {currentData.summary.L1_nonzero}
                        </a>{" "}
                        (out of {currentData.summary.L1_total})
                      </TableCell>
                      <TableCell>
                        <a href="#" style={{ color: "#1976d2" }}>
                          {currentData.summary.L2_nonzero}
                        </a>{" "}
                        (out of {currentData.summary.L2_total})
                      </TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell>
                        <strong>Mean ({activeFeature})</strong>
                      </TableCell>
                      <TableCell>{currentData.summary.L1_mean}</TableCell>
                      <TableCell>{currentData.summary.L2_mean}</TableCell>
                    </TableRow>
                    <TableRow>
                      <TableCell>
                        <strong>Median ({activeFeature})</strong>
                      </TableCell>
                      <TableCell>{currentData.summary.L1_median}</TableCell>
                      <TableCell>{currentData.summary.L2_median}</TableCell>
                    </TableRow>
                  </TableBody>
                </Table>
              </TableContainer>

              <Typography
                variant="subtitle2"
                sx={{
                  bgcolor: "#f5f5f5",
                  p: 1,
                  borderTop: "1px solid #ccc",
                  borderLeft: "1px solid #ccc",
                  borderRight: "1px solid #ccc",
                }}
              >
                ⚗ Statistical testing results: six p-values of the three
                statistical tests
              </Typography>
              <Box sx={{ border: "1px solid #ccc", borderTop: 0, p: 2 }}>
                <ul
                  style={{ fontSize: "0.85rem", color: "#555", marginTop: 0 }}
                >
                  <li>
                    Analyzed quantitative feature (QF):{" "}
                    <strong>{activeFeature}</strong>
                  </li>
                  <li>
                    QF (L1) &gt; QF (L2): The QF in L1 is larger than the QF in
                    L2
                  </li>
                  <li>
                    The p-values which are less than the{" "}
                    <strong>
                      {method} correction: p-value cutoff = {cutoffStr}
                    </strong>{" "}
                    are highlighted with the{" "}
                    <span style={{ backgroundColor: "yellow" }}>
                      yellow background
                    </span>
                    .
                  </li>
                </ul>

                <TableContainer
                  component={Paper}
                  variant="outlined"
                  sx={{ width: "80%", margin: "0 auto", borderRadius: 0 }}
                >
                  <Table size="small">
                    <TableHead sx={{ bgcolor: "#fafafa" }}>
                      <TableRow>
                        <TableCell></TableCell>
                        <TableCell>
                          <strong>t-test</strong>
                        </TableCell>
                        <TableCell>
                          <strong>U test</strong>
                        </TableCell>
                        <TableCell>
                          <strong>KS test</strong>
                        </TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {currentData.tests.map((row, idx) => {
                        const tVal = getPValue(row, "t");
                        const uVal = getPValue(row, "u");
                        const ksVal = getPValue(row, "ks");
                        return (
                          <TableRow key={idx}>
                            <TableCell>
                              <strong>{row.direction}</strong>
                            </TableCell>
                            <TableCell
                              sx={{
                                bgcolor:
                                  tVal < cutoffValue ? "#fff59d" : "inherit",
                              }}
                            >
                              {formatPValue(tVal)}
                            </TableCell>
                            <TableCell
                              sx={{
                                bgcolor:
                                  uVal < cutoffValue ? "#fff59d" : "inherit",
                              }}
                            >
                              {formatPValue(uVal)}
                            </TableCell>
                            <TableCell
                              sx={{
                                bgcolor:
                                  ksVal < cutoffValue ? "#fff59d" : "inherit",
                              }}
                            >
                              {formatPValue(ksVal)}
                            </TableCell>
                          </TableRow>
                        );
                      })}
                    </TableBody>
                  </Table>
                </TableContainer>
              </Box>
            </Box>
          )}

          {/* ⭐ 圖表視圖 (補回 Plotly 程式碼) */}
          {viewMode === "figure" && currentPlot && (
            <Grid container spacing={4}>
              <Grid item xs={12} md={6}>
                <Box
                  sx={{
                    width: "100%",
                    height: 500,
                    border: "1px solid #eee",
                    p: 1,
                  }}
                >
                  <Plot
                    data={[
                      {
                        y: currentPlot.L1,
                        type: "box",
                        name: "L1",
                        boxpoints: "all",
                        jitter: 0.3,
                      },
                      {
                        y: currentPlot.L2,
                        type: "box",
                        name: "L2",
                        boxpoints: "all",
                        jitter: 0.3,
                      },
                    ]}
                    layout={{
                      title: "Box Plot",
                      autosize: true,
                      margin: { l: 50, r: 20, t: 50, b: 50 },
                    }}
                    useResizeHandler={true}
                    style={{ width: "100%", height: "100%" }}
                  />
                </Box>
              </Grid>
              <Grid item xs={12} md={6}>
                <Box
                  sx={{
                    width: "100%",
                    height: 500,
                    border: "1px solid #eee",
                    p: 1,
                  }}
                >
                  <Plot
                    data={[
                      {
                        ...generateCdfData(currentPlot.L1),
                        mode: "lines",
                        name: "L1",
                      },
                      {
                        ...generateCdfData(currentPlot.L2),
                        mode: "lines",
                        name: "L2",
                      },
                    ]}
                    layout={{
                      title: "CDF Plot",
                      xaxis: { title: activeFeature },
                      yaxis: { title: "Prob (X ≤ x)" },
                      autosize: true,
                      margin: { l: 50, r: 20, t: 50, b: 50 },
                    }}
                    useResizeHandler={true}
                    style={{ width: "100%", height: "100%" }}
                  />
                </Box>
              </Grid>
            </Grid>
          )}
        </Box>
      </Box>
    );
  };

  return (
    <Box
      sx={{
        p: 4,
        maxWidth: "1200px",
        margin: "0 auto",
        fontFamily: "sans-serif",
      }}
    >
      <Typography
        variant="h5"
        color="#1976d2"
        gutterBottom
        sx={{ display: "flex", alignItems: "center", mb: 3 }}
      >
        <DescriptionIcon sx={{ mr: 1, fontSize: 30 }} /> Compare Quantitative
        Features Between Two Yeast Gene Lists
      </Typography>
      {/* ⭐ Step 1 與 Step 2 並排的 Grid 結構 */}
      <Grid container spacing={3} sx={{ mb: 3 }}>
        {/* Step 1 區塊 (左半邊) */}
        <Grid item xs={12} md={6}>
          <Paper
            variant="outlined"
            sx={{ height: "100%", bgcolor: "#f8f9fa", borderRadius: 1 }}
          >
            <Box
              sx={{
                p: 2,
                borderBottom: "1px solid #ddd",
                display: "flex",
                justifyContent: "space-between",
                alignItems: "center",
              }}
            >
              <Typography variant="subtitle1" fontWeight="bold">
                Step1. Input Two Gene Lists (L1 vs. L2)
              </Typography>
              {/* 也可以在這裡放一個全局的 Load All Sample 按鈕，如果你想要的話 */}
            </Box>
            <Box sx={{ p: 2 }}>
              <Grid container spacing={2}>
                <Grid item xs={6}>
                  <Box
                    sx={{
                      border: "1px solid #ddd",
                      p: 1.5,
                      bgcolor: "white",
                      borderRadius: 1,
                      height: "100%",
                    }}
                  >
                    <Typography variant="body2" fontWeight="bold" mb={1}>
                      Input First Gene List (L1) <br />
                      {/* ⭐ 改成可點擊的 span */}
                      <span
                        style={{
                          color: "#1976d2",
                          fontSize: "0.8rem",
                          fontWeight: "normal",
                          cursor: "pointer",
                          textDecoration: "underline",
                        }}
                        onClick={handleLoadSample1}
                      >
                        Load Sample
                      </span>
                    </Typography>
                    <TextField
                      multiline
                      rows={10}
                      fullWidth
                      value={list1Text}
                      onChange={(e) => setList1Text(e.target.value)}
                      sx={{ bgcolor: "white" }}
                    />
                  </Box>
                </Grid>
                <Grid item xs={6}>
                  <Box
                    sx={{
                      border: "1px solid #ddd",
                      p: 1.5,
                      bgcolor: "white",
                      borderRadius: 1,
                      height: "100%",
                    }}
                  >
                    <Typography variant="body2" fontWeight="bold" mb={1}>
                      Input Second Gene List (L2) <br />
                      {/* ⭐ 改成可點擊的 span */}
                      <span
                        style={{
                          color: "#1976d2",
                          fontSize: "0.8rem",
                          fontWeight: "normal",
                          cursor: "pointer",
                          textDecoration: "underline",
                        }}
                        onClick={handleLoadSample2}
                      >
                        Load Sample
                      </span>
                    </Typography>
                    <TextField
                      multiline
                      rows={10}
                      fullWidth
                      value={list2Text}
                      onChange={(e) => setList2Text(e.target.value)}
                      sx={{ bgcolor: "white" }}
                    />
                  </Box>
                </Grid>
              </Grid>
            </Box>
          </Paper>
        </Grid>

        {/* Step 2 區塊 (右半邊) */}
        <Grid item xs={12} md={6}>
          <Paper
            variant="outlined"
            sx={{ height: "100%", bgcolor: "#f8f9fa", borderRadius: 1 }}
          >
            <Box sx={{ p: 2, borderBottom: "1px solid #ddd" }}>
              <Typography variant="subtitle1" fontWeight="bold">
                Step2. Select the Quantitative Feature to be Analyzed
              </Typography>
            </Box>
            <Box
              sx={{
                p: 2,
                bgcolor: "white",
                m: 2,
                border: "1px solid #ddd",
                borderRadius: 1,
                height: "350px",
                overflowY: "auto",
              }}
            >
              {FEATURE_GROUPS.map((group) => (
                <Box key={group.category} sx={{ mb: 2 }}>
                  <Typography
                    variant="body2"
                    color="#1976d2"
                    fontWeight="bold"
                    sx={{ mb: 1 }}
                  >
                    ⊕ {group.features.length} {group.category}
                  </Typography>
                  <FormGroup sx={{ ml: 3 }}>
                    {group.features.map((feat) => (
                      <FormControlLabel
                        key={feat}
                        control={
                          <Checkbox
                            size="small"
                            name={feat}
                            checked={selectedFeatures[feat] || false}
                            onChange={handleFeatureChange}
                          />
                        }
                        label={<Typography variant="body2">{feat}</Typography>}
                        sx={{ m: 0, height: 28 }}
                      />
                    ))}
                  </FormGroup>
                </Box>
              ))}
            </Box>
          </Paper>
        </Grid>
      </Grid>{" "}
      {/* 結束 Step 1 & 2 並排區塊 */}
      {/* ⭐ Step 3 置於下方 */}
      <Paper
        variant="outlined"
        sx={{ bgcolor: "#f8f9fa", borderRadius: 1, mb: 3 }}
      >
        <Box sx={{ p: 2, borderBottom: "1px solid #ddd" }}>
          <Typography variant="subtitle1" fontWeight="bold">
            Step3. Specify P-value Cutoff for Multiple Hypothesis Testing
          </Typography>
        </Box>
        <Box sx={{ p: 2 }}>
          <RadioGroup
            value={method}
            onChange={(e) => setMethod(e.target.value)}
          >
            <FormControlLabel
              value="Bonferroni"
              control={<Radio size="small" />}
              label={
                <Box display="flex" alignItems="center">
                  <Typography variant="body2">
                    <strong>Bonferroni correction:</strong> p-value cutoff
                    ={" "}
                  </Typography>
                  <TextField
                    size="small"
                    variant="outlined"
                    value={method === "Bonferroni" ? cutoffStr : "10^-3"}
                    onChange={(e) => setCutoffStr(e.target.value)}
                    disabled={method !== "Bonferroni"}
                    sx={{
                      width: 90,
                      ml: 1,
                      "& input": { p: 0.5, textAlign: "center" },
                    }}
                  />
                </Box>
              }
              sx={{ mb: 1 }}
            />
            <FormControlLabel
              value="FDR"
              control={<Radio size="small" />}
              label={
                <Box display="flex" alignItems="center">
                  <Typography variant="body2">
                    FDR (False Discovery Rate): p-value cutoff ={" "}
                  </Typography>
                  <TextField
                    size="small"
                    variant="outlined"
                    value={method === "FDR" ? cutoffStr : "10^-3"}
                    onChange={(e) => setCutoffStr(e.target.value)}
                    disabled={method !== "FDR"}
                    sx={{
                      width: 90,
                      ml: 1,
                      "& input": { p: 0.5, textAlign: "center" },
                    }}
                  />
                </Box>
              }
              sx={{ mb: 1 }}
            />
            <FormControlLabel
              value="None"
              control={<Radio size="small" />}
              label={
                <Box display="flex" alignItems="center">
                  <Typography variant="body2">
                    No correction: p-value cutoff ={" "}
                  </Typography>
                  <TextField
                    size="small"
                    variant="outlined"
                    value={method === "None" ? cutoffStr : "10^-3"}
                    onChange={(e) => setCutoffStr(e.target.value)}
                    disabled={method !== "None"}
                    sx={{
                      width: 90,
                      ml: 1,
                      "& input": { p: 0.5, textAlign: "center" },
                    }}
                  />
                </Box>
              }
            />
          </RadioGroup>
        </Box>
      </Paper>
      {/* 按鈕置中 */}
      <Box sx={{ textAlign: "center" }}>
        <Button
          variant="outlined"
          color="warning"
          sx={{ mr: 2, textTransform: "none", fontWeight: "bold" }}
          onClick={() => {
            setList1Text("");
            setList2Text("");
            setSelectedFeatures({});
          }}
        >
          Reset ⟲
        </Button>
        <Button
          variant="contained"
          color="primary"
          onClick={handleAnalyze}
          disabled={loading}
          sx={{ textTransform: "none", minWidth: 120, fontWeight: "bold" }}
        >
          {loading ? (
            <CircularProgress size={24} color="inherit" />
          ) : (
            "Submit ⮞"
          )}
        </Button>
      </Box>
      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
      {renderResultSection()}
    </Box>
  );
}
