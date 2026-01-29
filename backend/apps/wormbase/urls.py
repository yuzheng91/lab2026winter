from django.urls import path
from . import views

app_name = "wormbase"

urlpatterns = [
    path('api/analyze/', views.analyze_worm, name='analyze_worm'),
]