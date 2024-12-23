from collections import namedtuple

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
        "saved",
        "savepath",
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
    cbar_shrink=1,
    cbar_aspect=25,
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
)

default_distribution_state = dict(
  visualization="histogram",
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
    showTitles=True,
    subtitle="Histogram",
    title="Histogram",
    xtitle="Percent Identity",
    ytitle="Frequency",
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
    showTitles=True,
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
    showTitles=True,
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
    saved=False,
    savepath="",
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
    return DocState(
        id=id,
        view=view,
        dataView=dataView,
        filename=filename,
        filetype=filetype,
        filemtime=filemtime,
        basename=basename,
        saved=saved,
        savepath=savepath,
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

    # state = default_state

    # def get_state():
    #     return state

    # def set_state(**kwargs):
    #     nonlocal state

    #     state = state._replace(**kwargs)
    #     on_state_updated()

    # def reset_state():
    #     nonlocal state
    #     state = default_state
    #     on_state_updated()

    # def on_state_updated():
    #     if on_update:
    #         on_update(state)

    # return get_state, set_state, reset_state
