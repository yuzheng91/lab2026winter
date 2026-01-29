import os
import json
import math
import numpy as np
import pandas as pd
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# 設定檔案路徑 (請根據你的專案結構調整)
APP_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(APP_DIR, "data")
ID_FILE = os.path.join(DATA_DIR, "c_elegans.PRJNA13758.WS298.geneOtherIDs.txt")
GAF_FILE = os.path.join(DATA_DIR, "gene_association.WS298.wb.c_elegans")

# === 全域變數 (伺服器啟動時載入，避免每次 Request 都重讀) ===
# 格式: { 'F': { 'GO:001': {genes}, 'GO:002': {genes} }, 'P': ..., 'C': ... }
_GO_DATABASE = {"F": {}, "P": {}, "C": {}} 
_LIVE_GENE_IDS = set() # 所有的 Live Gene WBGeneID
_SYMBOL_MAP = {}       # Symbol/Name -> WBGeneID

def load_data():
    """ 預先載入並處理 WormBase 資料 """
    global _LIVE_GENE_IDS, _SYMBOL_MAP, _GO_DATABASE
    
    if _LIVE_GENE_IDS: return # 已經載入過就跳過

    print("Loading WormBase Data...")

    # 1. 處理 Gene IDs 對照表 (過濾 Dead genes)
    # 假設檔案格式為 CSV/TSV，第一欄是 ID，第二欄是 Status (Live/Dead)
    # 這裡依照常見 WormBase 格式假設 header
    try:
        # 讀取 Gene IDs (根據你的檔案格式可能需要調整 sep 或 header)
        df_ids = pd.read_csv(ID_FILE, sep='\t', header=None, names=['GeneID', 'Status', 'SeqName', 'Name', 'Other'], encoding='utf-8')
        
        # 只留 Live
        live_df = df_ids[df_ids['Status'] == 'Live']
        _LIVE_GENE_IDS = set(live_df['GeneID'].values)

        # 建立名稱映射 (Name -> ID, SeqName -> ID)
        # 轉成大寫以利比對
        for _, row in live_df.iterrows():
            gid = row['GeneID']
            if pd.notna(row['Name']): _SYMBOL_MAP[str(row['Name']).upper()] = gid
            if pd.notna(row['SeqName']): _SYMBOL_MAP[str(row['SeqName']).upper()] = gid
            _SYMBOL_MAP[gid.upper()] = gid # 自己對自己也要 mapping
            
    except Exception as e:
        print(f"Error loading ID file: {e}")

    # 2. 處理 GAF 檔案 (Gene Association File)
    try:
        # GAF 是標準格式，通常用 tab 分隔，註解以 ! 開頭
        # 我們需要 Col 1 (DB_Object_ID), Col 4 (GO ID), Col 8 (Aspect) -> 0-based index: 1, 4, 8
        # 注意: Python read_csv 0-based index: 
        # Col 2: DB_Object_ID (WBGene...)
        # Col 5: GO ID (GO:...)
        # Col 9: Aspect (F, P, C)
        
        gaf_cols = [1, 4, 8] # 讀取這三欄即可
        df_gaf = pd.read_csv(GAF_FILE, sep='\t', comment='!', header=None, usecols=gaf_cols, names=['GeneID', 'GOID', 'Aspect'])

        # 只保留 Live Genes
        df_gaf = df_gaf[df_gaf['GeneID'].isin(_LIVE_GENE_IDS)]

        # 依照 Aspect 分組建立 Database
        # 結構: Aspect -> GO_ID -> Set of GeneIDs
        for aspect, group in df_gaf.groupby('Aspect'):
            # Group by GOID inside this aspect
            go_map = group.groupby('GOID')['GeneID'].apply(set).to_dict()
            _GO_DATABASE[aspect] = go_map

    except Exception as e:
        print(f"Error loading GAF file: {e}")
    
    print("WormBase Data Loaded.")

# 確保啟動時執行一次 (在 Django 中可以放在 AppConfig，這裡簡單呼叫)
# 注意：在開發模式下存檔自動重啟會觸發這個
if not _LIVE_GENE_IDS:
    load_data()


def _perform_enrichment(input_genes_set, aspect_char, method, threshold):
    """ 對單一 Aspect 執行分析 """
    
    go_map = _GO_DATABASE.get(aspect_char, {})
    
    # 背景母體數 (N): 該 Aspect 下所有被註解到的 Live Genes 聯集總數
    # 或者用所有 Live Genes 總數? 通常用 "有該 Aspect 註解的基因總數" 比較準
    background_genes_in_aspect = set().union(*go_map.values())
    F_total = len(background_genes_in_aspect) # N
    
    # 篩選 Input genes: 必須在背景母體中
    valid_input_genes = input_genes_set.intersection(background_genes_in_aspect)
    S_input = len(valid_input_genes) # n (Draw size)

    results = []

    if S_input == 0:
        return {"enriched": [], "depleted": []}

    for go_id, gene_set in go_map.items():
        # G: 背景中擁有該 GO term 的基因數 (Success in population)
        G_count = len(gene_set)
        
        # O: Input 中擁有該 GO term 的基因數 (Success in draw)
        overlap_genes = valid_input_genes.intersection(gene_set)
        O_count = len(overlap_genes)
        
        # 效能優化: 如果 Input 裡完全沒這個 GO term，Enriched p-value 必定是 1，可以跳過 Enriched 計算
        # 但 Depleted 還是要算
        
        # Hypergeometric Test
        # Enriched: P(X >= O) = sf(O-1)
        p_enriched = hypergeom.sf(O_count - 1, F_total, G_count, S_input)
        
        # Depleted: P(X <= O) = cdf(O)
        p_depleted = hypergeom.cdf(O_count, F_total, G_count, S_input)
        
        # 計算 Fold Change (Observed / Expected)
        expected = (G_count / F_total) * S_input
        raw_ratio = (O_count / expected) if expected > 0 else 0
        if raw_ratio > 0:
            fold_change = math.log2(raw_ratio)
        else:
            fold_change = 0

        results.append({
            "GO ID": go_id,
            "Gene Count": G_count, # 背景總數
            "Input Hit": O_count,  # input 中到的數
            "Expected": expected,
            "Fold Change": fold_change,
            "P-value Enriched": p_enriched,
            "P-value Depleted": p_depleted
        })

    # Convert to DataFrame for easier correction
    df = pd.DataFrame(results)
    if df.empty:
         return {"enriched": [], "depleted": []}

    # Multiple Testing Correction (Separately for Enriched and Depleted)
    # 1. Enriched
    p_vals_enr = df["P-value Enriched"].values
    if method == 'FDR':
        _, adj_p, _, _ = multipletests(p_vals_enr, method='fdr_bh')
        col_name = "FDR"
    elif method == 'Bonferroni':
        _, adj_p, _, _ = multipletests(p_vals_enr, method='bonferroni')
        col_name = "Bonferroni"
    else:
        adj_p = p_vals_enr
        col_name = "P-value" # No correction
    
    df["Adj P-value Enriched"] = adj_p

    # 2. Depleted
    p_vals_dep = df["P-value Depleted"].values
    if method == 'FDR':
        _, adj_p_dep, _, _ = multipletests(p_vals_dep, method='fdr_bh')
    elif method == 'Bonferroni':
        _, adj_p_dep, _, _ = multipletests(p_vals_dep, method='bonferroni')
    else:
        adj_p_dep = p_vals_dep
    
    df["Adj P-value Depleted"] = adj_p_dep

    # Filtering & Sorting
    thresh_val = float(threshold) if threshold != "None" else 1.0
    
    # Split into two tables
    df_enr_final = df[df["Adj P-value Enriched"] <= thresh_val].copy()
    df_enr_final = df_enr_final.sort_values("Adj P-value Enriched").head(500) # Limit rows
    
    df_dep_final = df[df["Adj P-value Depleted"] <= thresh_val].copy()
    df_dep_final = df_dep_final.sort_values("Adj P-value Depleted").head(500)

    # Helper to format output
    def format_records(dframe, p_col):
        out = []
        for _, row in dframe.iterrows():
            out.append({
                "GO ID": row["GO ID"],
                "Expected": f"{row['Gene Count']}/{F_total} ({row['Gene Count']/F_total:.2%})", # 這裡簡化顯示
                "Observed": f"{row['Input Hit']}/{S_input} ({row['Input Hit']/S_input:.2%})",
                "Fold Change": row["Fold Change"],
                "P-value": row[p_col], # 這是校正後的
                "Raw P-value": row["P-value Enriched" if "Enriched" in p_col else "P-value Depleted"]
            })
        return out

    return {
        "enriched": format_records(df_enr_final, "Adj P-value Enriched"),
        "depleted": format_records(df_dep_final, "Adj P-value Depleted")
    }

@csrf_exempt
def analyze_worm(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=405)

    try:
        payload = json.loads(request.body)
        text = payload.get("text", "")
        method = payload.get("method", "None")
        threshold = payload.get("threshold", "0.05")
        
        # 1. Parse Input & Map to IDs
        input_list = [line.strip().upper() for line in text.splitlines() if line.strip()]
        mapped_ids = set()
        for item in input_list:
            if item in _SYMBOL_MAP:
                mapped_ids.add(_SYMBOL_MAP[item])
        
        if not mapped_ids:
             return JsonResponse({"error": "No valid Live C. elegans genes found."}, status=400)

        # 2. Run Analysis for F, P, C
        response_data = {
            "meta": {
                "input_total": len(input_list),
                "mapped_live": len(mapped_ids)
            },
            "F": _perform_enrichment(mapped_ids, "F", method, threshold),
            "P": _perform_enrichment(mapped_ids, "P", method, threshold),
            "C": _perform_enrichment(mapped_ids, "C", method, threshold)
        }
        
        return JsonResponse(response_data)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)