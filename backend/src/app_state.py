# Usage
# get_state, set_state, reset_state = create_app_state()

# Update state
# set_state(filename="example.txt")
# set_state(progress=50)
# set_state(debug=True)

# Access state
# current_state = get_state()
# print("Current State:", current_state)

from collections import namedtuple

AppState = namedtuple(
    "AppState",
    [
        "view",
        "filename",
        "filetype",
        "basename",
        "progress",
        "stage",
        "debug",
        "tempdir_path",  # TODO: rename tempdir_path, it's no longer an accurate name
        "sequences_count",
        "alignment_output_path",
        "export_path",
        "pair_progress",
        "pair_count",
        "estimated_time",
        "validation_error_id",
        "compute_stats",
        "platform",
    ],
)


def create_app_state(
    view="runner",
    filename="",
    filetype="",
    basename="",
    progress=0,
    stage="",
    debug=False,
    on_update=None,
    tempdir_path="",
    sequences_count=0,
    alignment_output_path="",
    export_path="",
    pair_progress=0,
    pair_count=0,
    estimated_time=None,
    validation_error_id=None,
    compute_stats=None,
    platform=None,
):
    default_state = AppState(
        view=view,
        filename=filename,
        filetype=filetype,
        basename=basename,
        progress=progress,
        stage=stage,
        debug=debug,
        tempdir_path=tempdir_path,
        sequences_count=sequences_count,
        alignment_output_path=alignment_output_path,
        export_path=export_path,
        pair_progress=pair_progress,
        pair_count=pair_count,
        estimated_time=estimated_time,
        validation_error_id=validation_error_id,
        compute_stats=compute_stats,
        platform=platform,
    )

    state = default_state

    def get_state():
        return state

    def set_state(**kwargs):
        nonlocal state

        if state.debug == True:
            print(kwargs)

        state = state._replace(**kwargs)
        on_state_updated()

    def reset_state():
        nonlocal state
        state = default_state
        on_state_updated()

    def on_state_updated():
        if on_update:
            on_update(state)

    return get_state, set_state, reset_state
