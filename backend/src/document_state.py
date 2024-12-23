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
        "progress",
        "stage",
        "tempdir_path",
        "sequences_count",
        "pair_progress",
        "pair_count",
        "estimated_time",
        "validation_error_id",
        "compute_stats",
    ],
)


def create_document_state(
    id,
    view="runner",
    dataView="heatmap",
    filename="",
    filetype="",
    filemtime=None,
    basename="",
    progress=0,
    stage="",
    tempdir_path="",
    sequences_count=0,
    pair_progress=0,
    pair_count=0,
    estimated_time=None,
    validation_error_id=None,
    compute_stats=None,
):
    return DocState(
        id=id,
        view=view,
        dataView=dataView,
        filename=filename,
        filetype=filetype,
        filemtime=filemtime,
        basename=basename,
        progress=progress,
        stage=stage,
        tempdir_path=tempdir_path,
        sequences_count=sequences_count,
        pair_progress=pair_progress,
        pair_count=pair_count,
        estimated_time=estimated_time,
        validation_error_id=validation_error_id,
        compute_stats=compute_stats,
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
