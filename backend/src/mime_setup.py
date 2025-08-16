import mimetypes


def register_mimetypes():
    mimetypes.add_type("text/fasta", ".fasta")
    mimetypes.add_type("text/fasta", ".fas")
    mimetypes.add_type("text/fasta", ".faa")
    mimetypes.add_type("text/fasta", ".fnt")
    mimetypes.add_type("text/fasta", ".fa")
    mimetypes.add_type("application/vnd.sdt", ".sdt")
