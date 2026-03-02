import os, re
import json
import math
import numpy as np
import pandas as pd
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# ==========================================
# 1. 檔案路徑與全域變數設定
# ==========================================
APP_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(APP_DIR, "data")
ID_FILE = os.path.join(DATA_DIR, "c_elegans.PRJNA13758.WS298.geneOtherIDs.txt")
GAF_FILE = os.path.join(DATA_DIR, "gene_association.WS298.wb.c_elegans")
# 加入新的 Feature Table 路徑
FEATURE_TABLE_FILE = os.path.join(DATA_DIR, "d12_fgg_55feature_table.csv")

_GO_DATABASE = {"F": {}, "P": {}, "C": {}} 
_FEATURE_GROUP_DATABASE = {} # 存放 55 個 feature 的資料庫
_LIVE_GENE_IDS = set() 
_SYMBOL_MAP = {}

_EVIDENCE_DF = None 
_ID_TO_SYMBOL = {}

# ==========================================
# 2. 系統初始化與資料載入
# ==========================================
def load_data():
    global _LIVE_GENE_IDS, _SYMBOL_MAP, _GO_DATABASE, _FEATURE_GROUP_DATABASE
    global _EVIDENCE_DF, _ID_TO_SYMBOL # ★ 宣告全域變數

    if _LIVE_GENE_IDS: return 
    print("Loading Data...")

    # --- 載入 WormBase ID ---
    try:
        df_ids = pd.read_csv(ID_FILE, sep='\t', header=None, names=['GeneID', 'Status', 'SeqName', 'Name', 'Other'], encoding='utf-8')
        live_df = df_ids[df_ids['Status'] == 'Live']
        _LIVE_GENE_IDS = set(live_df['GeneID'].values)
        for _, row in live_df.iterrows():
            gid = row['GeneID']
            # ★ 新增：建立 ID 對應 Symbol 的字典，給 Gene Groups 顯示用
            if pd.notna(row['Name']): 
                _SYMBOL_MAP[str(row['Name']).upper()] = gid
                _ID_TO_SYMBOL[gid] = str(row['Name'])
            elif pd.notna(row['SeqName']):
                _SYMBOL_MAP[str(row['SeqName']).upper()] = gid
                _ID_TO_SYMBOL[gid] = str(row['SeqName'])
            else:
                _ID_TO_SYMBOL[gid] = gid
            _SYMBOL_MAP[gid.upper()] = gid 
    except Exception as e:
        print(f"Error loading ID file: {e}")

    # --- 載入 GAF FPC 資料 ---
    try:
        # ★ 修改：原本只讀 [1, 4, 8]，現在多讀 [2, 5, 6] 拿 Symbol, Reference, Evidence
        gaf_cols = [1, 2, 4, 5, 6, 8]
        df_gaf = pd.read_csv(GAF_FILE, sep='\t', comment='!', header=None, usecols=gaf_cols, 
                             names=['GeneID', 'Symbol', 'GOID', 'Reference', 'Evidence', 'Aspect'])
        df_gaf = df_gaf[df_gaf['GeneID'].isin(_LIVE_GENE_IDS)]
        
        # ★ 新增：將完整的 DataFrame 存起來給 API 查詢用
        _EVIDENCE_DF = df_gaf 

        for aspect, group in df_gaf.groupby('Aspect'):
            go_map = group.groupby('GOID')['GeneID'].apply(set).to_dict()
            _GO_DATABASE[aspect] = go_map
    except Exception as e:
        print(f"Error loading GAF file: {e}")
        
    # ... (保留你的 55 Feature 資料載入邏輯) ...

    print("All Data Loaded.")

if not _LIVE_GENE_IDS:
    load_data()


# ==========================================
# 3. Pipeline A: 原版 FPC 分析
# ==========================================
def _perform_fpc_enrichment(input_genes_set, aspect_char):
    """ 原本的 FPC 分析邏輯，保持不動 """
    go_map = _GO_DATABASE.get(aspect_char, {})
    F_total = len(_LIVE_GENE_IDS)
    
    valid_input_genes = input_genes_set.intersection(_LIVE_GENE_IDS)
    S_input = len(valid_input_genes)

    if S_input == 0:
         return {"enriched": [], "depleted": []}

    results = []
    for go_id, gene_set in go_map.items():
        G_count = len(gene_set)
        overlap_genes = valid_input_genes.intersection(gene_set)
        O_count = len(overlap_genes)
        
        p_enriched = hypergeom.sf(O_count - 1, F_total, G_count, S_input)
        p_depleted = hypergeom.cdf(O_count, F_total, G_count, S_input)
        
        expected = (G_count / F_total) * S_input
        raw_ratio = (O_count / expected) if expected > 0 else 0
        fold_change = math.log2(raw_ratio) if raw_ratio > 0 else 0

        results.append({
            "GO ID": go_id, "Gene Count": G_count, "Input Hit": O_count,
            "Expected": expected, "Fold Change": fold_change,
            "raw_p_enr": p_enriched, "raw_p_dep": p_depleted
        })

    df = pd.DataFrame(results)
    if df.empty:
         return {"enriched": [], "depleted": []}

    for dir_suffix in ["enr", "dep"]:
        p_vals = df[f"raw_p_{dir_suffix}"].values
        _, fdr_vals, _, _ = multipletests(p_vals, method='fdr_bh')
        _, bon_vals, _, _ = multipletests(p_vals, method='bonferroni')
        df[f"fdr_p_{dir_suffix}"] = fdr_vals
        df[f"bon_p_{dir_suffix}"] = bon_vals

    def format_records(dframe, dir_suffix):
        out = []
        for _, row in dframe.iterrows():
            out.append({
                "GO ID": row["GO ID"], "Gene Count": row["Gene Count"],
                "Input Hit": row["Input Hit"], "Expected": row["Expected"],
                "Fold Change": row["Fold Change"], "TotalInAspect": F_total, "TotalInput": S_input,
                "p_raw": row[f"raw_p_{dir_suffix}"], "p_fdr": row[f"fdr_p_{dir_suffix}"], "p_bon": row[f"bon_p_{dir_suffix}"]
            })
        return out

    df_enr = df.sort_values("raw_p_enr")
    df_dep = df.sort_values("raw_p_dep")
    return {"enriched": format_records(df_enr, "enr"), "depleted": format_records(df_dep, "dep")}


# ==========================================
# 4. Pipeline B: 新版 55 Feature 分析 (來自 test.py)
# ==========================================
def _perform_gene_group_enrichment(input_genes_set):
    """ 實作 test.py 的邏輯 (不鎖死背景母體數) """
    N = len(_LIVE_GENE_IDS)
    valid_input_genes = input_genes_set.intersection(_LIVE_GENE_IDS)
    n = len(valid_input_genes)

    if n == 0:
        return {"enriched": [], "depleted": []}

    results = []

    for feature_name, feature_genes in _FEATURE_GROUP_DATABASE.items():
        M = len(feature_genes)
        overlap_genes = valid_input_genes & feature_genes
        k = len(overlap_genes)
        
        expected_ratio = M / N if N > 0 else 0
        observed_ratio = k / n if n > 0 else 0
        fold_change = (observed_ratio / expected_ratio) if expected_ratio > 0 else 0
        
        # 字串格式化
        expected_count = (n * M) / N
        expected_ratio_str = f"{round(expected_count, 2)}/{n} ({round(expected_ratio * 100, 2)}%)"
        observed_ratio_str = f"{k}/{n} ({round(observed_ratio * 100, 2)}%)"
        
        p_enriched = hypergeom.sf(k - 1, N, M, n)
        p_depleted = hypergeom.cdf(k, N, M, n)

        results.append({
            "Feature Name": feature_name,
            "Observed Ratio": observed_ratio_str,
            "Expected Ratio": expected_ratio_str,
            "Fold Change": fold_change,
            "raw_p_enr": p_enriched,
            "raw_p_dep": p_depleted
        })
        
    df = pd.DataFrame(results)
    if df.empty:
        return {"enriched": [], "depleted": []}

    # 計算 FDR 與 Bonferroni
    for dir_suffix in ["enr", "dep"]:
        p_vals = df[f"raw_p_{dir_suffix}"].values
        _, fdr_vals, _, _ = multipletests(p_vals, method='fdr_bh')
        _, bon_vals, _, _ = multipletests(p_vals, method='bonferroni')
        df[f"fdr_p_{dir_suffix}"] = fdr_vals
        df[f"bon_p_{dir_suffix}"] = bon_vals

    def format_records(dframe, dir_suffix):
        out = []
        for _, row in dframe.iterrows():
            fc = row["Fold Change"]
            log2_fc = math.log2(fc) if fc > 0 else -1000
            
            pval = row[f"raw_p_{dir_suffix}"]
            log10_pval = -math.log10(pval) if pval > 0 else 1000
            
            fdr_pval = row[f"fdr_p_{dir_suffix}"]
            log10_fdr = -math.log10(fdr_pval) if fdr_pval > 0 else 1000

            out.append({
                "Feature Name": row["Feature Name"],
                "Observed Ratio": row["Observed Ratio"],
                "Expected Ratio": row["Expected Ratio"],
                "Fold Change": fc,
                "P-value": pval,
                "FDR P-value": fdr_pval,
                "Bonferroni P-value": row[f"bon_p_{dir_suffix}"],
                "log2(Fold Change)": log2_fc,
                "-log10(P-value)": log10_pval,
                "-log10(FDR P-value)": log10_fdr
            })
        return out

    df_enr = df.sort_values("raw_p_enr")
    df_dep = df.sort_values("raw_p_dep")

    return {
        "enriched": format_records(df_enr, "enr"),
        "depleted": format_records(df_dep, "dep")
    }


# ==========================================
# 5. API 進入點 (分發器 Router)
# ==========================================
@csrf_exempt
def analyze_worm(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=405)

    try:
        payload = json.loads(request.body)
        text = payload.get("text", "")
        
        # 前端傳來的開關 (預設跑原本的，保護舊系統)
        run_fpc = payload.get("run_fpc", True) 
        run_gene_group = payload.get("run_gene_group", False)
        
        input_list = [line.strip().upper() for line in text.splitlines() if line.strip()]
        mapped_ids = set()
        for item in input_list:
            if item in _SYMBOL_MAP:
                mapped_ids.add(_SYMBOL_MAP[item])
        
        if not mapped_ids:
             return JsonResponse({"error": "No valid Live C. elegans genes found."}, status=400)

        # 準備基礎回傳資料
        response_data = {
            "meta": {
                "input_total": len(input_list),
                "mapped_live": len(mapped_ids),
                "total_background": len(_LIVE_GENE_IDS)
            }
        }
        
        # 根據前端勾選，掛載分析結果
        if run_fpc:
            response_data["fpc_results"] = {
                "F": _perform_fpc_enrichment(mapped_ids, "F"),
                "P": _perform_fpc_enrichment(mapped_ids, "P"),
                "C": _perform_fpc_enrichment(mapped_ids, "C")
            }
            
        if run_gene_group:
            response_data["gene_group_results"] = _perform_gene_group_enrichment(mapped_ids)
            
        return JsonResponse(response_data)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)
    
@csrf_exempt
def fetch_term_evidence(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=405)

    global _LIVE_GENE_IDS, _EVIDENCE_DF
    if not _LIVE_GENE_IDS or _EVIDENCE_DF is None:
        print("[Debug] ⚠️ _EVIDENCE_DF 為空，正在重新載入資料...")
        load_data()

    try:
        payload = json.loads(request.body)
        term_id = payload.get("term_id", "").strip()
        text = payload.get("genes", "")
        
        print(f"\n========== 🔍 準備查詢明細 ==========")
        print(f"1. 前端要求的 GO ID: [{term_id}]")
        
        raw_genes = re.split(r'[\n\r,\t\s]+', text)
        input_list = [g.strip().upper() for g in raw_genes if g.strip()]
        print(f"2. 前端傳來的基因數量 (去空白後): {len(input_list)}")
        
        mapped_ids = set()
        for item in input_list:
            if item in _SYMBOL_MAP:
                mapped_ids.add(_SYMBOL_MAP[item])
                
        print(f"3. 成功轉換成 WormBase ID 的數量: {len(mapped_ids)}")
        if len(mapped_ids) > 0:
            print(f"   (範例 ID: {list(mapped_ids)[:3]}...)")
        
        rows = []
        if term_id.startswith("GO:"):
            if _EVIDENCE_DF is not None:
                print(f"4. _EVIDENCE_DF 總筆數: {len(_EVIDENCE_DF)}")
                
                # 分步驟檢查，看是哪邊沒有交集
                term_match_count = (_EVIDENCE_DF['GOID'].str.strip() == term_id).sum()
                print(f"5. 資料庫中符合 [{term_id}] 的筆數: {term_match_count}")
                
                gene_match_count = _EVIDENCE_DF['GeneID'].str.strip().isin(mapped_ids).sum()
                print(f"6. 資料庫中符合輸入基因的筆數: {gene_match_count}")
                
                # 實際取交集
                mask = (_EVIDENCE_DF['GOID'].str.strip() == term_id) & (_EVIDENCE_DF['GeneID'].str.strip().isin(mapped_ids))
                res_df = _EVIDENCE_DF[mask].copy()
                print(f"7. 最終交集 (GO + Gene) 命中筆數 (去重前): {len(res_df)}")
                
                res_df = res_df.drop_duplicates(subset=['GeneID', 'GOID'])
                print(f"8. 最終回傳筆數 (去重後): {len(res_df)}")
                
                for _, row in res_df.iterrows():
                    symbol = row["Symbol"]
                    if pd.isna(symbol): symbol = row["GeneID"]
                    rows.append({
                        "WormBase_ID": str(row["GeneID"]).strip(),
                        "Symbol": str(symbol).strip(),
                        "GO_ID": str(row["GOID"]).strip(),
                        "Reference": str(row["Reference"]).strip() if pd.notna(row["Reference"]) else "-",
                        "Evidence_Code": str(row["Evidence"]).strip() if pd.notna(row["Evidence"]) else "-",
                    })
            else:
                print("[Debug] ❌ _EVIDENCE_DF 依然是 None！載入失敗！")
        
        print("======================================\n")
        return JsonResponse({"rows": rows})

    except Exception as e:
        import traceback
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)