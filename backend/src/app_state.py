# Usage
# get_state, set_state = create_app_state()

# Update state
# set_state(filename="example.txt")
# set_state(progress=50)
# set_state(debug=True)

# Access state
# current_state = get_state()
# print("Current State:", current_state)

from collections import namedtuple

AppState = namedtuple('AppState', [
    'view',
    'filename',
    'filetype',
    'basename',
    'progress',
    'debug',
    'tempdir_path', # TODO: rename tempdir_path, it's no longer an accurate name
    'sequences_count',
    'alignment_output_path',
    'export_path',
    'performance_profiles'
])

def create_app_state(
    view="runner",
    filename="",
    filetype="",
    basename="",
    progress=0,
    debug=False,
    on_update=None,
    tempdir_path="",
    sequences_count=0,
    alignment_output_path='',
    export_path='',
    performance_profiles=dict(),
):
    default_state = AppState(
        view=view,
        filename=filename,
        filetype=filetype,
        basename=basename,
        progress=progress,
        debug=debug,
        tempdir_path=tempdir_path,
        sequences_count=sequences_count,
        alignment_output_path=alignment_output_path,
        export_path=export_path,
        performance_profiles=performance_profiles
    )

    state = default_state

    def get_state():
        return state


    def set_state(**kwargs):
        nonlocal state
        state = state._replace(**kwargs)
        on_state_updated()

    def reset_state():
        nonlocal state
        state = default_state
        on_state_updated()

    def on_state_updated():
        if state.debug == True:
            print("State updated:", state)

        if on_update:
            on_update(state)

    return get_state, set_state, reset_state
