import os
import math
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# ==========================================
# 1. 設定與路徑
# ==========================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
ID_FILE = os.path.join(DATA_DIR, "c_elegans.PRJNA13758.WS298.geneOtherIDs.txt")
GAF_FILE = os.path.join(DATA_DIR, "gene_association.WS298.wb.c_elegans")

# 全域資料庫
_GO_DATABASE = {"F": {}, "P": {}, "C": {}} 
_LIVE_GENE_IDS = set() 
_SYMBOL_MAP = {} 

# 強制設定背景母體數 (依照老師的結果)
FORCE_TOTAL_GENES = 49164 

# ==========================================
# 2. 資料載入
# ==========================================
def load_data():
    global _LIVE_GENE_IDS, _SYMBOL_MAP, _GO_DATABASE
    
    if _LIVE_GENE_IDS: return
    print(f"[System] Loading data...")

    # A. 載入 ID 對照 (用於 Symbol 轉換)
    try:
        if os.path.exists(ID_FILE):
            df_ids = pd.read_csv(ID_FILE, sep='\t', header=None, 
                                 names=['GeneID', 'Status', 'SeqName', 'Name', 'Other'], 
                                 encoding='utf-8')
            live_df = df_ids[df_ids['Status'] == 'Live']
            _LIVE_GENE_IDS = set(live_df['GeneID'].values)

            for _, row in live_df.iterrows():
                gid = row['GeneID']
                if pd.notna(row['Name']): _SYMBOL_MAP[str(row['Name']).upper()] = gid
                if pd.notna(row['SeqName']): _SYMBOL_MAP[str(row['SeqName']).upper()] = gid
                _SYMBOL_MAP[gid.upper()] = gid 
        else:
            print("[Warning] ID file not found. Mapping might fail.")
    except Exception as e:
        print(f"[Error] ID Load Failed: {e}")

    # B. 載入 GAF
    try:
        if os.path.exists(GAF_FILE):
            gaf_cols = [1, 4, 8]
            df_gaf = pd.read_csv(GAF_FILE, sep='\t', comment='!', header=None, 
                                 usecols=gaf_cols, names=['GeneID', 'GOID', 'Aspect'])
            
            # 確保只用 Live Genes (如果有載入的話)
            if _LIVE_GENE_IDS:
                df_gaf = df_gaf[df_gaf['GeneID'].isin(_LIVE_GENE_IDS)]

            for aspect, group in df_gaf.groupby('Aspect'):
                go_map = group.groupby('GOID')['GeneID'].apply(set).to_dict()
                _GO_DATABASE[aspect] = go_map
    except Exception as e:
        print(f"[Error] GAF Load Failed: {e}")

# ==========================================
# 3. 核心運算 (仿照老師格式)
# ==========================================
def _perform_enrichment_teacher_style(input_genes_set, aspect_char):
    
    go_map = _GO_DATABASE.get(aspect_char, {})
    
    # 1. 設定母體 N (老師使用的是 49164)
    F_total = FORCE_TOTAL_GENES
    
    # 2. 篩選 Input (n)
    # 這裡假設使用者輸入的 32 個基因都是 Valid 的 (老師範例中 Observed Ratio 分母為 32)
    # 如果 input_genes_set 包含非 Live 基因，這裡會過濾掉
    valid_input_genes = input_genes_set.intersection(_LIVE_GENE_IDS) if _LIVE_GENE_IDS else input_genes_set
    S_input = len(valid_input_genes)

    if S_input == 0:
        return None, None

    results = []
    
    # 遍歷所有 GO Terms
    for go_id, gene_set in go_map.items():
        G_count = len(gene_set) # K
        
        overlap_genes = valid_input_genes.intersection(gene_set)
        O_count = len(overlap_genes) # k (Observed)
        
        # 計算 P-values
        # Enriched: P(X >= k)
        p_enriched = hypergeom.sf(O_count - 1, F_total, G_count, S_input)
        
        # Depleted: P(X <= k)
        p_depleted = hypergeom.cdf(O_count, F_total, G_count, S_input)
        
        # 計算數值
        # Expected Ratio (Bg Ratio)
        exp_ratio_str = f"{G_count}/{F_total} ({G_count/F_total:.4%})"
        
        # Observed Ratio
        obs_ratio_str = f"{O_count}/{S_input} ({O_count/S_input:.4%})" # 老師格式顯示較多位數? 範例為 3.125%
        # 修正百分比顯示格式以匹配範例 "0.002%"
        def fmt_pct(num, den):
            pct = (num / den) * 100
            # 簡單邏輯：如果是整數就顯示整數，否則顯示小數
            return f"{num}/{den} ({pct:.4g}%)"

        exp_ratio_str = fmt_pct(G_count, F_total)
        obs_ratio_str = fmt_pct(O_count, S_input)

        # Fold Change (Raw)
        expected = (G_count / F_total) * S_input
        fc_raw = (O_count / expected) if expected > 0 else 0
        
        # log2 FC (老師用 -1000.0 代表 0)
        fc_log2 = math.log2(fc_raw) if fc_raw > 0 else -1000.0

        results.append({
            "Domain Name": go_id,
            "Expected Ratio": exp_ratio_str,
            "Observed Ratio": obs_ratio_str,
            "Enrichment Fold Change": fc_raw,
            "log2(Enrichment Fold Change)": fc_log2,
            "P-value Enriched": p_enriched,
            "P-value Depleted": p_depleted
        })

    df = pd.DataFrame(results)
    if df.empty: return pd.DataFrame(), pd.DataFrame()

    # --- 製作 Enriched 表 ---
    df_enr = df.copy()
    # P-value 欄位
    df_enr["P-value"] = df_enr["P-value Enriched"]
    # FDR / Bonferroni
    _, fdr_enr, _, _ = multipletests(df_enr["P-value"], method='fdr_bh')
    _, bon_enr, _, _ = multipletests(df_enr["P-value"], method='bonferroni')
    df_enr["FDR enriched P-value"] = fdr_enr
    df_enr["Bonferroni enriched P-value"] = bon_enr
    
    # -log10 欄位
    df_enr["-log10(P-value)"] = -np.log10(df_enr["P-value"] + 1e-300) # 避免 log(0)
    df_enr["-log10(FDR P-value)"] = -np.log10(df_enr["FDR enriched P-value"] + 1e-300)
    df_enr["-log10(Bonferroni P-value)"] = -np.log10(df_enr["Bonferroni enriched P-value"] + 1e-300)
    
    # 修正 -0.0 顯示 (老師檔案有 -0.0)
    # Python 預設 float 處理即可
    
    # --- 製作 Depleted 表 ---
    df_dep = df.copy()
    df_dep["P-value"] = df_dep["P-value Depleted"]
    _, fdr_dep, _, _ = multipletests(df_dep["P-value"], method='fdr_bh')
    _, bon_dep, _, _ = multipletests(df_dep["P-value"], method='bonferroni')
    df_dep["FDR depleted P-value"] = fdr_dep
    df_dep["Bonferroni depleted P-value"] = bon_dep

    df_dep["-log10(P-value)"] = -np.log10(df_dep["P-value"] + 1e-300)
    df_dep["-log10(FDR P-value)"] = -np.log10(df_dep["FDR depleted P-value"] + 1e-300)
    df_dep["-log10(Bonferroni P-value)"] = -np.log10(df_dep["Bonferroni depleted P-value"] + 1e-300)

    # 選擇欄位並排序 (依照老師格式)
    cols_enr = [
        "Domain Name", "Expected Ratio", "Observed Ratio", "Enrichment Fold Change",
        "P-value", "FDR enriched P-value", "Bonferroni enriched P-value",
        "log2(Enrichment Fold Change)", "-log10(P-value)", "-log10(FDR P-value)", "-log10(Bonferroni P-value)"
    ]
    cols_dep = [
        "Domain Name", "Expected Ratio", "Observed Ratio", "Enrichment Fold Change",
        "P-value", "FDR depleted P-value", "Bonferroni depleted P-value",
        "log2(Enrichment Fold Change)", "-log10(P-value)", "-log10(FDR P-value)", "-log10(Bonferroni P-value)"
    ]
    
    # 老師的 Depleted 檔案似乎是依照 P-value 排序? (snippet 0.999 在前?) 
    # 其實通常是依照 Domain Name 或 P-value。這裡我們預設依照 P-value 排序。
    # 註：如果想要完全一樣的順序，可能需要知道老師的排序邏輯。這裡我們用 Domain Name 排序以方便對照。
    return df_enr[cols_enr].sort_values("Domain Name"), df_dep[cols_dep].sort_values("Domain Name")

# ==========================================
# 4. 主程式
# ==========================================
def main(input_genes_list):
    load_data()
    
    # 映射 Input
    mapped_ids = set()
    for item in input_genes_list:
        clean = item.strip().upper()
        if clean in _SYMBOL_MAP: mapped_ids.add(_SYMBOL_MAP[clean])
    
    print(f"Mapped Genes: {len(mapped_ids)}")
    if len(mapped_ids) == 0: return

    OUTPUT_DIR = "output_results"
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

    aspect_names = {"F": "Molecular_Function", "P": "Biological_Process", "C": "Cellular_Component"}
    
    for code, name in aspect_names.items():
        print(f"Processing {name} ({code})...")
        df_enr, df_dep = _perform_enrichment_teacher_style(mapped_ids, code)
        
        if df_enr is not None:
            # 存檔 - Enriched
            enr_name = f"GO_Term_{code}_enriched_result.csv"
            df_enr.to_csv(os.path.join(OUTPUT_DIR, enr_name), index=False, encoding='utf-8-sig')
            
            # 存檔 - Depleted
            dep_name = f"GO_Term_{code}_depleted_result.csv"
            df_dep.to_csv(os.path.join(OUTPUT_DIR, dep_name), index=False, encoding='utf-8-sig')
            
            print(f"  -> Saved {enr_name} and {dep_name}")

    print("Done.")

if __name__ == "__main__":
    MY_SAMPLE_GENES = [
        "WBGene00002324", "WBGene00003164", "WBGene00004343", "WBGene00007012",
        "WBGene00007013", "WBGene00007017", "WBGene00007018", "WBGene00007026",
        "WBGene00020951", "WBGene00001078", "WBGene00007189", "WBGene00010720",
        "WBGene00017983", "WBGene00018175", "WBGene00020375", "WBGene00020820",
        "WBGene00003582", "WBGene00004879", "WBGene00004880", "WBGene00004881",
        "WBGene00004882", "WBGene00004883", "WBGene00004884", "WBGene00004885",
        "WBGene00010551", "WBGene00012343", "WBGene00012904", "WBGene00013188",
        "WBGene00015943", "WBGene00018156", "WBGene00021365", "WBGene00021929"
    ]
    main(MY_SAMPLE_GENES)