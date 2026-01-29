# enrichment_analysis/views.py
import os, json
import numpy as np
import pandas as pd
import math
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


APP_DIR = os.path.dirname(__file__)
FILE_PATH = os.path.join(APP_DIR, "data", "domain_test.xlsx")

_domain_df = pd.read_excel(FILE_PATH)

_domains = _domain_df["Domains"].astype(str).to_numpy()
_C = _domain_df["number"].astype(int).to_numpy()
_id_sets = _domain_df["ID"].fillna("").astype(str).map(lambda s: set(s.split(","))).tolist()
_F = 6705


def _compute_overlap_counts(input_set: set[str]) -> np.ndarray:
    return np.fromiter((len(input_set & s) for s in _id_sets), dtype=int, count=len(_id_sets))


def _enrichment_pvals(T: np.ndarray, S: int, G: np.ndarray, F: int):
    # 超幾何分佈：X ~ Hypergeom(M=F, n=G, N=S)
    # enriched: P(X >= T) = sf(T-1)
    # depleted: P(X <= T) = cdf(T)
    p_greater = hypergeom.sf(T - 1, F, G, S)
    p_less = hypergeom.cdf(T, F, G, S)
    return p_greater, p_less

def index(request):
    return render(request, "enrichment_analysis/index.html")


@csrf_exempt
def analyze(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=405)

    try:
        payload = json.loads(request.body.decode("utf-8"))
    except Exception:
        return JsonResponse({"error": "Invalid JSON"}, status=400)

    text = (payload.get("text") or "").strip()
    if not text:
        return JsonResponse({"error": "text is required"}, status=400)
    filter_method = payload.get("method", "None")
    filter_threshold = payload.get("threshold", "0.05")

    # 1) parse input list
    input_list = [line.strip() for line in text.splitlines() if line.strip()]
    input_set = set(input_list)
    B_total = float(len(input_list))  # input genes count
    if B_total == 0:
        return JsonResponse({"error": "No valid inputs"}, status=400)

    # 2) overlap counts (A_arr / T) and vectorized pvals
    T = _compute_overlap_counts(input_set)         # overlap counts per domain
    A_arr = np.array(T, dtype=float)
    C_arr = np.array(_C, dtype=float)
    D_total = float(_F)                             # 6705
    N = len(_domains)

    pvalue_greater_list, pvalue_less_list = _enrichment_pvals(T=T, S=int(B_total), G=_C, F=_F)

    # 3) corrections (keep exactly the objects you asked for)
    CUT_OFF = 0.01
    P_value_corr_enriched_FDR = multipletests(pvalue_greater_list, alpha=CUT_OFF, method="fdr_bh")
    P_value_corr_depleted_FDR = multipletests(pvalue_less_list, alpha=CUT_OFF, method="fdr_bh")
    P_value_corr_enriched_Bonferroni = multipletests(pvalue_greater_list, alpha=CUT_OFF, method="bonferroni")
    P_value_corr_depleted_Bonferroni = multipletests(pvalue_less_list, alpha=CUT_OFF, method="bonferroni")

    # 4) ratios + fold change (same format you specified)
    expected_pct = (C_arr / D_total) * 100.0
    observed_pct = (A_arr / B_total) * 100.0

    df_enriched_result = pd.DataFrame({
        "Domain Name": _domains.astype(str),
        "Expected Ratio": [f"{int(c)}/{int(D_total)} ({p:.4f}%)" for c, p in zip(C_arr, expected_pct)],
        "Observed Ratio": [f"{int(a)}/{int(B_total)} ({p:.4f}%)" for a, p in zip(A_arr, observed_pct)],
        "Fold Change": [
            math.log2(float(o)/float(e)) if (float(e) > 0 and float(o) > 0) 
            else float("-inf") if float(o) == 0 and float(e) > 0
            else float("inf") 
            for e, o in zip(expected_pct, observed_pct)
        ],
        "P-value": pvalue_greater_list[:N],
        "FDR enriched P-value": P_value_corr_enriched_FDR[1][:N],
        "Bonferroni enriched P-value": P_value_corr_enriched_Bonferroni[1][:N],
    })

    df_depleted_result = pd.DataFrame({
        "Domain Name": _domains.astype(str),
        "Expected Ratio": [f"{int(c)}/{int(D_total)} ({p:.4f}%)" for c, p in zip(C_arr, expected_pct)],
        "Observed Ratio": [f"{int(a)}/{int(B_total)} ({p:.4f}%)" for a, p in zip(A_arr, observed_pct)],
        "Fold Change": [
            math.log2(float(o)/float(e)) if (float(e) > 0 and float(o) > 0) 
            else float("-inf") if float(o) == 0 and float(e) > 0
            else float("inf") 
            for e, o in zip(expected_pct, observed_pct)
        ],
        "P-value": pvalue_less_list[:N],
        "FDR depleted P-value": P_value_corr_depleted_FDR[1][:N],
        "Bonferroni depleted P-value": P_value_corr_depleted_Bonferroni[1][:N],
    })

    print(df_enriched_result.head())
    print(len(df_enriched_result))
    print(df_depleted_result.head())
    print(len(df_depleted_result))

    # 5) Normalize column names & Backend Filtering
    # 統一欄位名稱，讓前端 DataGrid 可以共用 columns 設定
    rename_map_enriched = {
        "FDR enriched P-value": "FDR",
        "Bonferroni enriched P-value": "Bonferroni"
    }
    rename_map_depleted = {
        "FDR depleted P-value": "FDR",
        "Bonferroni depleted P-value": "Bonferroni"
    }
    df_enriched_result.rename(columns=rename_map_enriched, inplace=True)
    df_depleted_result.rename(columns=rename_map_depleted, inplace=True)

    # df_out = df.sort_values("P-value", ascending=True).head(500)    # (optional but safe) sort by enriched p-value so frontend looks sane

    df_enriched_result = df_enriched_result.replace(
    [np.inf, -np.inf], None).where(df_enriched_result.notna(), None)
    df_depleted_result = df_depleted_result.replace(
    [np.inf, -np.inf], None).where(df_depleted_result.notna(), None)

    return JsonResponse({
        "meta": {
            "input_count": int(B_total),
            "total_domains": int(N),
            "enriched_count": len(df_enriched_result),
            "depleted_count": len(df_depleted_result),
        },
        "enriched": df_enriched_result.to_dict(orient="records"),
        "depleted": df_depleted_result.to_dict(orient="records"),
    })