from django.urls import path
from . import views

app_name = "yeast"

urlpatterns = [
    path("api/analyze/", views.analyze, name="analyze"),
]