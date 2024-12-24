from collections import namedtuple
from document_state import DocState, create_document_state, load_doc_settings

AppState = namedtuple(
    "AppState",
    [
        "debug",
        "platform",
        "documents"
    ],
)


def create_app_state(
    debug=False,
    on_update=None,
    platform=None,
    documents=[]
):
    default_state = AppState(
        debug=debug,
        platform=platform,
        documents=documents
    )

    state = default_state

    def get_state():
        return state

    def set_state(updater=None, **kwargs):
        nonlocal state
        current_state = state
        if callable(updater):
           state = state._replace(**updater(current_state)._asdict())
        else:
           state = state._replace(**kwargs)
           on_state_updated()

    # def set_state(skip_callbacks=False, **kwargs):
    #     nonlocal state
    #     state = state._replace(**kwargs)

    #     if skip_callbacks == False:
    #         on_state_updated()

    def reset_state():
        nonlocal state
        state = default_state
        on_state_updated()

    def on_state_updated():
        if on_update:
            on_update(state)

    def new_document(id: str, **kwargs):
        doc = get_document(id)
        if doc:
            doc_settings = load_doc_settings(kwargs["tempdir_path"])
            if doc_settings:
                return update_document(id, **kwargs, **doc_settings)
            else:
                return update_document(id, **kwargs)

        doc_state = create_document_state(id=id, **kwargs)
        set_state(documents=state.documents + [doc_state])

    def get_document(id: str) -> DocState | None:
        return next((doc for doc in state.documents if doc.id == id), None)

    def update_document(id: str, **updates):
        updated_documents = [
            doc._replace(**updates) if doc.id == id else doc
            for doc in state.documents
        ]
        set_state(documents=updated_documents)

    def remove_document(id: str):
        updated_documents = [doc for doc in state.documents if doc.id != id]
        set_state(documents=updated_documents)

    def remove_empty_documents():
        updated_documents = [doc for doc in state.documents if doc.filename != ""]
        set_state(documents=updated_documents)

    return get_state, set_state, reset_state, new_document, get_document, update_document, remove_document, remove_empty_documents
