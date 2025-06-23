# import re
# from Bio.SeqUtils import gc_fraction
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord


# def preprocess(raw_seq_dict):
#     seq_dict = {}
#     for key, record in raw_seq_dict.items():
#         sequence = str(record.seq).upper()
#         sequence = re.sub(r"[^A-Z]", "", sequence)
#         key= re.sub(r"[, -]", "_", key)
#         key = key[:50] ## hardcap for id length
#         seq_dict[key] = str(sequence)
#     return seq_dict

# def grab_stats(seq_dict):
#     seq_stats = {}
#     for key, record_str in seq_dict.items():
#         gcCount = round(gc_fraction(record_str), 2) * 100
#         genLen = len(record_str)
#         seq_stats.setdefault(key, [])
#         seq_stats[key].append(gcCount)
#         seq_stats[key].append(genLen)
#     return seq_stats

# def residue_check(seq):
#     return bool(re.search(r"[EFILPQZ]", seq))

# def preprocess_fasta_to_fasta(input_fasta_path, output_fasta_path):
#     processed_records = []
#     for record in SeqIO.parse(input_fasta_path, "fasta"):
#         # Process ID
#         processed_id = re.sub(r"[, -]", "_", record.id)
#         processed_id = processed_id[:50] # hardcap for id length

#         # Process sequence
#         sequence_str = str(record.seq).upper()
#         sequence_str = re.sub(r"[^A-Z]", "", sequence_str)

#         # Create new SeqRecord
#         new_record = SeqRecord(Seq(sequence_str), id=processed_id, description="") # Clear description
#         processed_records.append(new_record)

#     SeqIO.write(processed_records, output_fasta_path, "fasta")
#     return output_fasta_path
