from django.urls import path
from . import views

app_name = "wormbase"

urlpatterns = [
    # 這裡只需要定義這個 app 自己的路徑就好
    path('api/analyze/', views.analyze_worm, name='analyze_worm'),
    path('api/fetch_term_evidence/', views.fetch_term_evidence, name='fetch_term_evidence'),
    
    # ⚠️ 請確保這份檔案裡面絕對沒有 include('wormbase.urls') 這一行！
]