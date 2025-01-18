from collections import namedtuple
import json
import os

DocState = namedtuple(
    "DocState",
    [
        "id",
        "view",
        "dataView",
        "filename",
        "filetype",
        "filemtime",
        "basename",
        "modified",
        "progress",
        "stage",
        "tempdir_path",
        "sequences_count",
        "pair_progress",
        "pair_count",
        "estimated_time",
        "validation_error_id",
        "compute_stats",
        "heatmap",
        "distribution"
    ],
)

default_heatmap_state = dict(
    colorScaleKey="Portland",
    reverse=False,
    vmax=100,
    vmin=65,
    cellspace=1,
    annotation=False,
    annotation_font_size=10,
    annotation_rounding=0,
    annotation_alpha="0",
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
    axis_labels=False,
    axlabel_xrotation=270,
    axlabel_xfontsize=12,
    axlabel_yrotation=360,
    axlabel_yfontsize=12,
    cutoff_1=95,
    cutoff_2=75
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
  visualization="histogram", # TODO: all of these need to be removed
  dataSet="scores",
  histogram=dict(
    **visualization_defaults,
    binColor="hsl(195, 53%, 79%)",
    binSize=1,
    histOutlineWidth=1,
    histlineColor="hsl(9, 100%, 64%)",
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
  ),
  raincloud=dict(
    **visualization_defaults,
    bandwidth=8,
    jitter=0.5,
    markerColor="hsl(9, 100%, 64%)",
    markerSize=7,
    plotOrientation="horizontal",
    pointPos=-1.5,
    points="all",
    showAxisLines=True,
    showPoints=True,
    showZeroLine=False,
    fillColor="hsl(195, 53%, 79%)",
    editable=False,
    side="positive",
    showMeanline=True,
    makeEditable=True,
    dticks=5,
    subtitle="Raincloud Plot",
    title="Raincloud Plot",
    xtitle="Percent Identity",
    ytitle="Genome",
  ),
  violin=dict(
    **visualization_defaults,
    bandwidth=5,
    boxWidth=0.5,
    boxfillColor="hsl(195, 53%, 79%)",
    boxlineColor="hsl(9, 100%, 64%)",
    boxlineWidth=3,
    fillColor="hsl(195, 53%, 79%)",
    jitter=0.5,
    markerColor="hsl(9, 100%, 64%)",
    markerSize=7,
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

def create_document_state(
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
    heatmap=default_heatmap_state,
    distribution=default_distribution_state
):
    if filetype == "application/vnd.sdt" and tempdir_path:
        doc_settings = load_doc_settings(tempdir_path)
        if doc_settings:
            dataView = doc_settings["dataView"]
            heatmap = doc_settings["heatmap"]
            distribution = doc_settings["distribution"]

    return DocState(
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
        heatmap=heatmap,
        distribution=distribution
    )

def get_doc_setting_path(tempdir_path: str):
    return os.path.join(tempdir_path, "document.json")

def load_doc_settings(dir_path: str):
    path = get_doc_setting_path(dir_path)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None

def save_doc_settings(doc_state: DocState):
    path = get_doc_setting_path(doc_state.tempdir_path)
    with open(path, "w") as f:
        settings = {
            "dataView": doc_state.dataView,
            "heatmap": doc_state.heatmap,
            "distribution": doc_state.distribution
        }
        json.dump(settings, f, indent=2)
