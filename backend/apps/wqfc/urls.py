# backend/apps/wqfc/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('run_pipeline/', views.api_run_pipeline, name='wqfc_run_pipeline'),
    path('gene_detail/', views.api_gene_detail, name='wqfc_gene_detail'),
]