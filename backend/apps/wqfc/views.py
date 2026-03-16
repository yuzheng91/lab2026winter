import os
import json
import pandas as pd
import numpy as np

from django.conf import settings
from django.http import JsonResponse
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests
from django.views.decorators.csrf import csrf_exempt

@csrf_exempt
def api_run_pipeline(request):
    if request.method != "POST":
        return JsonResponse({"error": "Method not allowed. Please use POST."}, status=405)

    try:
        # 1. 解析前端傳來的 JSON
        body = json.loads(request.body)
        L1_raw = body.get("L1", "")
        L2_raw = body.get("L2", "")
        features = body.get("features", [])
        
        # 新增：接收前端傳來的檢定方法與閾值，雖然後端目前全算，但未來若需依賴此參數可預留
        correction_method = body.get("method", "Bonferroni") 
        p_cutoff = float(body.get("cutoff", 0.001))

        if not features:
            return JsonResponse({"status": "error", "message": "No features selected."}, status=400)

        L1 = [g.strip() for g in L1_raw.splitlines() if g.strip()]
        L2 = [g.strip() for g in L2_raw.splitlines() if g.strip()]

        # 2. 讀取核心資料檔
        excel_path = os.path.join(
            settings.BASE_DIR,
            "apps", "wqfc", "data", "Genome_features.xlsx"
        )
        if not os.path.exists(excel_path):
             return JsonResponse({"status": "error", "message": "Data file not found on server."}, status=500)
             
        df = pd.read_excel(excel_path).set_index("Feature name")

        # 建立別名映射字典 (這部分建議未來如果資料量大，可以在系統啟動時 load 進記憶體快取)
        alias_map = {}
        for idx, row in df.iterrows():
            if pd.notna(row.get("Alias")):
                for a in str(row["Alias"]).split("|"):
                    alias_map[a.strip()] = idx
            if pd.notna(row.get("Standard gene name")):
                alias_map[row["Standard gene name"]] = idx

        # 基因解析器
        def resolve_gene_list(gene_list):
            resolved, missing = [], []
            for g in gene_list:
                if g in df.index:
                    resolved.append(g)
                elif g in alias_map:
                    resolved.append(alias_map[g])
                else:
                    missing.append(g)
            return resolved, missing

        L1_resolved, missing_L1 = resolve_gene_list(L1)
        L2_resolved, missing_L2 = resolve_gene_list(L2)

        request.session["L1"] = L1
        request.session["L2"] = L2

        detail_dir = os.path.join(settings.BASE_DIR, "app", "static", "app", "detail")
        os.makedirs(detail_dir, exist_ok=True)

        plot_data, summaries, raw_tests_list = {}, {}, []

        # 3. 第一階段：計算基本統計量與原始檢定 p-value
        for feature in features:
            if feature not in df.columns:
                continue

            v1 = df.loc[df.index.intersection(L1_resolved), feature].dropna()
            v2 = df.loc[df.index.intersection(L2_resolved), feature].dropna()

            if len(v1) < 2 or len(v2) < 2:
                continue

            plot_data[feature] = {"L1": v1.tolist(), "L2": v2.tolist()}

            # 產生 Detail CSV (注意：高併發下這樣寫會有 Race Condition，建議檔名加上 session_id 或 timestamp)
            subset_L1 = df.loc[df.index.intersection(L1_resolved), ["Standard gene name", feature]].dropna()
            subset_L2 = df.loc[df.index.intersection(L2_resolved), ["Standard gene name", feature]].dropna()
            subset_L1.to_csv(os.path.join(detail_dir, f"L1_{feature}.csv"))
            subset_L2.to_csv(os.path.join(detail_dir, f"L2_{feature}.csv"))

            summaries[feature] = {
                "feature": feature,
                "L1_nonzero": int((v1 != 0).sum()),
                "L1_total": len(L1),
                "L2_nonzero": int((v2 != 0).sum()),
                "L2_total": len(L2),
                "L1_mean": f"{v1.mean():.3f}",     # 直接格式化為字串，確保前端顯示一致
                "L2_mean": f"{v2.mean():.3f}",
                "L1_median": f"{v1.median():.3f}",
                "L2_median": f"{v2.median():.3f}"
            }

            # 統計檢定
            t_stat, t_p_two = ttest_ind(v1, v2, equal_var=False)
            u_p_greater = mannwhitneyu(v1, v2, alternative="greater")[1]
            u_p_less = mannwhitneyu(v1, v2, alternative="less")[1]
            ks_p_greater = ks_2samp(v1, v2, alternative="less")[1]
            ks_p_less = ks_2samp(v1, v2, alternative="greater")[1]

            if v1.mean() > v2.mean():
                t_p_greater, t_p_less = t_p_two / 2, 1 - (t_p_two / 2)
            else:
                t_p_greater, t_p_less = 1 - (t_p_two / 2), t_p_two / 2

            raw_tests_list.extend([
                {"feature": feature, "direction": "QF (L1) > QF (L2)", "t_p": t_p_greater, "u_p": u_p_greater, "ks_p": ks_p_greater},
                {"feature": feature, "direction": "QF (L1) < QF (L2)", "t_p": t_p_less, "u_p": u_p_less, "ks_p": ks_p_less}
            ])

        results = {}

        # 4. 第二階段：執行全域多重檢驗校正
        if raw_tests_list:
            df_tests = pd.DataFrame(raw_tests_list)

            for test_type in ["t", "u", "ks"]:
                df_tests[f"{test_type}_fdr"] = multipletests(df_tests[f"{test_type}_p"], method="fdr_bh")[1]
                df_tests[f"{test_type}_bonf"] = multipletests(df_tests[f"{test_type}_p"], method="bonferroni")[1]

            # 5. 組裝 JSON
            for feature in summaries.keys():
                feat_tests = df_tests[df_tests["feature"] == feature]
                tests_output = []
                for _, row in feat_tests.iterrows():
                    tests_output.append({
                        "direction": row["direction"],
                        "t_test": float(row["t_p"]), "t_fdr": float(row["t_fdr"]), "t_bonf": float(row["t_bonf"]),
                        "u_test": float(row["u_p"]), "u_fdr": float(row["u_fdr"]), "u_bonf": float(row["u_bonf"]),
                        "ks_test": float(row["ks_p"]), "ks_fdr": float(row["ks_fdr"]), "ks_bonf": float(row["ks_bonf"])
                    })
                results[feature] = {"summary": summaries[feature], "tests": tests_output}

        summary_global = {
            "L1_total": len(L1),
            "L2_total": len(L2),
            "method": correction_method,
            "cutoff": p_cutoff
        }

        return JsonResponse({
            "status": "success",
            "summary_global": summary_global,
            "results": results,
            "plot_data": plot_data,
        })

    except Exception as e:
        import traceback
        traceback.print_exc() # 開發階段印出詳細錯誤
        return JsonResponse({"status": "error", "message": str(e)}, status=500)

# api_gene_detail 保持你的原樣即可


def api_gene_detail(request):
    if request.method != "GET":
        return JsonResponse({"error": "Method not allowed"}, status=405)
        
    feature = request.GET.get("feature")
    list_type = request.GET.get("list")

    file_path = os.path.join(
        settings.BASE_DIR,
        "apps", "wqfc", "data", "detail", f"{list_type}_{feature}.csv"
    )
    
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        return JsonResponse({"error": "Detail file not found"}, status=404)

    L1 = request.session.get("L1", [])
    L2 = request.session.get("L2", [])
    total = len(L1) if list_type == "L1" else len(L2)

    rows = []
    for _, row in df.iterrows():
        rows.append({
            "systematic": row.get("Feature name", ""),
            "standard": row.get("Standard gene name", ""),
            "value": row.get(feature, 0)
        })

    response = {
        "count": len(rows),
        "total": total,
        "feature": feature,
        "list": list_type,
        "data": rows
    }

    return JsonResponse(response)