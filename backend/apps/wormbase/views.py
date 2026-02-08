import os
import json
import math
import numpy as np
import pandas as pd
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# 設定檔案路徑
APP_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(APP_DIR, "data")
ID_FILE = os.path.join(DATA_DIR, "c_elegans.PRJNA13758.WS298.geneOtherIDs.txt")
GAF_FILE = os.path.join(DATA_DIR, "gene_association.WS298.wb.c_elegans")

# === 全域變數 ===
_GO_DATABASE = {"F": {}, "P": {}, "C": {}} 
_LIVE_GENE_IDS = set() 
_SYMBOL_MAP = {}

def load_data():
    """ 預先載入並處理 WormBase 資料 """
    global _LIVE_GENE_IDS, _SYMBOL_MAP, _GO_DATABASE
    if _LIVE_GENE_IDS: return 

    print("Loading WormBase Data...")
    try:
        df_ids = pd.read_csv(ID_FILE, sep='\t', header=None, names=['GeneID', 'Status', 'SeqName', 'Name', 'Other'], encoding='utf-8')
        live_df = df_ids[df_ids['Status'] == 'Live']
        _LIVE_GENE_IDS = set(live_df['GeneID'].values)
        for _, row in live_df.iterrows():
            gid = row['GeneID']
            if pd.notna(row['Name']): _SYMBOL_MAP[str(row['Name']).upper()] = gid
            if pd.notna(row['SeqName']): _SYMBOL_MAP[str(row['SeqName']).upper()] = gid
            _SYMBOL_MAP[gid.upper()] = gid 
    except Exception as e:
        print(f"Error loading ID file: {e}")

    try:
        gaf_cols = [1, 4, 8]
        df_gaf = pd.read_csv(GAF_FILE, sep='\t', comment='!', header=None, usecols=gaf_cols, names=['GeneID', 'GOID', 'Aspect'])
        df_gaf = df_gaf[df_gaf['GeneID'].isin(_LIVE_GENE_IDS)]
        for aspect, group in df_gaf.groupby('Aspect'):
            go_map = group.groupby('GOID')['GeneID'].apply(set).to_dict()
            _GO_DATABASE[aspect] = go_map
    except Exception as e:
        print(f"Error loading GAF file: {e}")
    print("WormBase Data Loaded.")

if not _LIVE_GENE_IDS:
    load_data()

def _perform_enrichment(input_genes_set, aspect_char):
    """ 
    對單一 Aspect 執行分析 
    (移除 method/threshold 參數，改為回傳所有數據)
    """
    go_map = _GO_DATABASE.get(aspect_char, {})
    F_total = len(_LIVE_GENE_IDS) # N
    
    valid_input_genes = input_genes_set.intersection(_LIVE_GENE_IDS)
    S_input = len(valid_input_genes) # n

    if S_input == 0:
        return {"enriched": [], "depleted": []}

    results = []

    # 1. 遍歷所有 GO ID (包含 Hit=0 的項目，照你要求保留)
    for go_id, gene_set in go_map.items():
        G_count = len(gene_set)
        
        overlap_genes = valid_input_genes.intersection(gene_set)
        O_count = len(overlap_genes)
        
        # Hypergeometric Test
        p_enriched = hypergeom.sf(O_count - 1, F_total, G_count, S_input)
        p_depleted = hypergeom.cdf(O_count, F_total, G_count, S_input)
        
        # Fold Change
        expected = (G_count / F_total) * S_input
        raw_ratio = (O_count / expected) if expected > 0 else 0
        fold_change = math.log2(raw_ratio) if raw_ratio > 0 else 0

        results.append({
            "GO ID": go_id,
            "Gene Count": G_count,
            "Input Hit": O_count,
            "Expected": expected,
            "Fold Change": fold_change,
            # 先存原始 P 值
            "raw_p_enr": p_enriched,
            "raw_p_dep": p_depleted
        })

    df = pd.DataFrame(results)
    if df.empty:
         return {"enriched": [], "depleted": []}

    # 2. 計算所有校正 P 值 (FDR & Bonferroni)
    
    # Enriched
    p_vals_enr = df["raw_p_enr"].values
    _, fdr_enr, _, _ = multipletests(p_vals_enr, method='fdr_bh')
    _, bon_enr, _, _ = multipletests(p_vals_enr, method='bonferroni')
    df["fdr_p_enr"] = fdr_enr
    df["bon_p_enr"] = bon_enr

    # Depleted
    p_vals_dep = df["raw_p_dep"].values
    _, fdr_dep, _, _ = multipletests(p_vals_dep, method='fdr_bh')
    _, bon_dep, _, _ = multipletests(p_vals_dep, method='bonferroni')
    df["fdr_p_dep"] = fdr_dep
    df["bon_p_dep"] = bon_dep

    # 3. 格式化回傳 (包含所有數據，不過濾)
    # 為了前端方便計算比例，把 TotalInAspect 和 TotalInput 放進去
    def format_records(dframe, dir_suffix):
        out = []
        for _, row in dframe.iterrows():
            out.append({
                "GO ID": row["GO ID"],
                "Gene Count": row["Gene Count"],
                "Input Hit": row["Input Hit"],
                "Expected": row["Expected"],
                "Fold Change": row["Fold Change"],
                "TotalInAspect": F_total,
                "TotalInput": S_input,
                
                # 回傳三種 P 值供前端篩選
                "p_raw": row[f"raw_p_{dir_suffix}"],
                "p_fdr": row[f"fdr_p_{dir_suffix}"],
                "p_bon": row[f"bon_p_{dir_suffix}"]
            })
        return out

    # 稍微排序 (依 Raw P)，方便前端預設顯示
    df_enr = df.sort_values("raw_p_enr")
    df_dep = df.sort_values("raw_p_dep")

    return {
        "enriched": format_records(df_enr, "enr"),
        "depleted": format_records(df_dep, "dep")
    }

@csrf_exempt
def analyze_worm(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=405)

    try:
        payload = json.loads(request.body)
        text = payload.get("text", "")
        # method/threshold 不再使用，因為後端全算
        
        input_list = [line.strip().upper() for line in text.splitlines() if line.strip()]
        mapped_ids = set()
        for item in input_list:
            if item in _SYMBOL_MAP:
                mapped_ids.add(_SYMBOL_MAP[item])
        
        if not mapped_ids:
             return JsonResponse({"error": "No valid Live C. elegans genes found."}, status=400)

        response_data = {
            "meta": {
                "input_total": len(input_list),
                "mapped_live": len(mapped_ids),
                "total_background": len(_LIVE_GENE_IDS)
            },
            "F": _perform_enrichment(mapped_ids, "F"),
            "P": _perform_enrichment(mapped_ids, "P"),
            "C": _perform_enrichment(mapped_ids, "C")
        }
        
        return JsonResponse(response_data)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return JsonResponse({"error": str(e)}, status=500)