"""Pydantic models for the data-extractor skill."""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field


class PlotType(str, Enum):
    SCATTER = "scatter"
    BAR = "bar"
    LINE = "line"
    BOX = "box"
    VIOLIN = "violin"
    HISTOGRAM = "histogram"
    HEATMAP = "heatmap"
    FOREST = "forest"
    KAPLAN_MEIER = "kaplan_meier"
    DOT_STRIP = "dot_strip"
    STACKED_BAR = "stacked_bar"
    FUNNEL = "funnel"
    ROC = "roc"
    VOLCANO = "volcano"
    WATERFALL = "waterfall"
    BLAND_ALTMAN = "bland_altman"
    PAIRED = "paired"
    BUBBLE = "bubble"
    AREA = "area"
    DOSE_RESPONSE = "dose_response"
    MANHATTAN = "manhattan"
    CORRELATION_MATRIX = "correlation_matrix"
    ERROR_BAR = "error_bar"
    TABLE = "table"
    OTHER = "other"
    NON_DATA = "non_data"


class Confidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class Figure(BaseModel):
    figure_id: str
    paper_id: str
    page_number: int
    image_index: int
    width: int
    height: int
    image_base64: str
    title: Optional[str] = None
    legend: Optional[str] = None
    plot_type: PlotType = PlotType.OTHER
    plot_type_confidence: Confidence = Confidence.MEDIUM


class DataSeries(BaseModel):
    name: str = "Series 1"
    x_values: list[float | str] = Field(default_factory=list)
    y_values: list[float] = Field(default_factory=list)
    error_bars_lower: list[float | None] = Field(default_factory=list)
    error_bars_upper: list[float | None] = Field(default_factory=list)


class ExtractedData(BaseModel):
    figure_id: str
    plot_type: PlotType
    title: Optional[str] = None
    x_label: Optional[str] = None
    y_label: Optional[str] = None
    x_unit: Optional[str] = None
    y_unit: Optional[str] = None
    x_min: Optional[float] = None
    x_max: Optional[float] = None
    y_min: Optional[float] = None
    y_max: Optional[float] = None
    x_scale: Optional[str] = None  # "linear" or "log"
    y_scale: Optional[str] = None  # "linear" or "log"
    series: list[DataSeries] = Field(default_factory=list)
    confidence: Confidence = Confidence.MEDIUM
    notes: Optional[str] = None
    legend_text: Optional[str] = None
    text_mentions: list[str] = Field(default_factory=list)
