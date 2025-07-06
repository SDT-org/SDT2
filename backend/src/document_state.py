from collections import namedtuple
import json
import os
from document_paths import build_document_paths

VERSION = 1

DocState = namedtuple(
    "DocState",
    [
        "version",
        "id",
        "view",
        "dataView",
        "filename",
        "filetype",
        "filemtime",
        "basename",
        "modified",
        "progress",  # TODO: remove
        "stage",  # TODO: remove
        "tempdir_path",
        "sequences_count",
        "pair_progress",
        "pair_count",
        "estimated_time",
        "validation_error_id",  # TODO: remove
        "compute_stats",
        "cluster_method",
        "heatmap",
        "clustermap",
        "distribution",
    ],
)

default_heatmap_state = dict(
    colorScaleKey="Portland",
    reverse=False,
    vmax=100,
    vmin=65,
    cellspace=1,
    annotation=True,
    annotation_rounding=2,
    showscale=True,
    titleFont="Sans Serif",
    showTitles=False,
    title="",
    subtitle="",
    xtitle="",
    ytitle="",
    cbar_shrink=5,
    cbar_aspect=2.5,
    cbar_pad=10,
    axis_labels=True,
    axlabel_xrotation=0,
    axlabel_fontsize=12,
    axlabel_yrotation=0,
    cutoff_1=95,
    cutoff_2=75,
)

default_clustermap_state = dict(
    threshold=70,
    method="average",
    annotation=True,
    titleFont="Sans Serif",
    showTitles=False,
    title="",
    axis_labels=True,
    axlabel_xrotation=0,
    axlabel_fontsize=12,
    axlabel_yrotation=0,
    cellspace=1,
)

visualization_defaults = dict(
    plotTitle="Distribution of Percent Identities",
    lineColor="hsl(9, 100%, 64%)",
    lineWidth=3,
    showAxisLabels=True,
    showGrid=True,
    showTickLabels=True,
    showTitles=True,
    titleFont="Sans Serif",
)

default_distribution_state = dict(
    visualization="distribution_histogram",
    dataSet="scores",
    histogram=dict(
        **visualization_defaults,
        binColor="hsl(195, 53%, 79%)",
        binSize=1,
        histOutlineWidth=0,
        histlineColor="hsl(0, 0%, 0%)",
        histnorm="probability",
        showHistogram=True,
        showLine=True,
        makeEditable=True,
        showAxisLines=True,
        showMeanline=True,
        dtickx=5,
        dticky=1,
        subtitle="Histogram",
        title="Histogram",
        xtitle="Percent Identity",
        ytitle="Frequency",
        plotOrientation="vertical",
        barGap=0.1,
    ),
    violin=dict(
        **visualization_defaults,
        bandwidth=5,
        boxWidth=0.05,
        boxfillColor="hsl(195, 53%, 79%)",
        boxlineColor="hsl(9, 100%, 64%)",
        boxlineWidth=3,
        fillColor="hsl(195, 53%, 79%)",
        jitter=0.5,
        markerColor="hsl(9, 100%, 64%)",
        markerSize=3,
        plotOrientation="vertical",
        pointOrientation="Violin",
        pointPos=0,
        points="all",
        showAxisLines=True,
        showBox=True,
        showMeanline=True,
        makeEditable=True,
        showPoints=True,
        showViolin=True,
        showZeroLine=False,
        whiskerWidth=0.2,
        title="Violin Plot",
        subtitle="Violin Plot",
        xtitle="",
        ytitle="",
    ),
)


def create_doc_state(
    id,
    view="runner",
    dataView="heatmap",
    filename="",
    filetype="",
    filemtime=None,
    basename="",
    modified=False,
    progress=0,
    stage="",
    tempdir_path="",
    sequences_count=0,
    pair_progress=0,
    pair_count=0,
    estimated_time=None,
    validation_error_id=None,
    compute_stats=None,
    cluster_method="average",
    heatmap=default_heatmap_state,
    clustermap=default_clustermap_state,
    distribution=default_distribution_state,
):
    if filetype == "application/vnd.sdt" and tempdir_path:
        settings = load_document_settings(tempdir_path)
        if settings:
            dataView = settings["dataView"]
            heatmap = settings["heatmap"]
            distribution = settings["distribution"]

    return DocState(
        version=VERSION,
        id=id,
        view=view,
        dataView=dataView,
        filename=filename,
        filetype=filetype,
        filemtime=filemtime,
        basename=basename,
        modified=modified,
        progress=progress,
        stage=stage,
        tempdir_path=tempdir_path,
        sequences_count=sequences_count,
        pair_progress=pair_progress,
        pair_count=pair_count,
        estimated_time=estimated_time,
        validation_error_id=validation_error_id,
        compute_stats=compute_stats,
        cluster_method=cluster_method,
        heatmap=heatmap,
        clustermap=clustermap,
        distribution=distribution,
    )


def load_document_settings(dir_path: str):
    doc_settings_path = build_document_paths(dir_path).settings
    if os.path.exists(doc_settings_path):
        with open(doc_settings_path) as f:
            return json.load(f)
    return None


def save_document_settings(doc_state: DocState):
    path = build_document_paths(doc_state.tempdir_path).settings
    with open(path, "w") as f:
        settings = {
            "version": VERSION,
            "dataView": doc_state.dataView,
            "heatmap": doc_state.heatmap,
            "clustermap": doc_state.clustermap,
            "distribution": doc_state.distribution,
        }
        json.dump(settings, f, indent=2)
