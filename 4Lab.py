import subprocess
from Bio import Entrez, SeqIO
from io import StringIO

def download_sequences(email, hpv_type, min_length=100):
    Entrez.email = email
    query = f"{hpv_type}[title] AND L1 [gene]" # Hpv type in L1 gene region
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1000)
    record = Entrez.read(handle)
    handle.close()

    if len(record['IdList']) > 0:
        sequence_ids = record['IdList']
        for sequence_id in sequence_ids:
            fasta_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            fasta_record = fasta_handle.read()
            fasta_handle.close()

            sequence = SeqIO.read(StringIO(fasta_record), "fasta")

            # Check sequence quality - if it's too short, skip it
            if len(sequence) < min_length:
                print(f"Skipping sequence {sequence.id} due to low quality")
                continue

            with open(f"HPV_{hpv_type}.fasta", "a") as outfile:
                SeqIO.write(sequence, outfile, "fasta")
    
    print(f"Downloaded sequences for HPV {hpv_type}")

def run_cd_hit_wsl(input_file, output_file, c=0.9):
    cmd = f"wsl cd-hit -i {input_file} -o {output_file} -c {c}" # 90% similarity threshold
    subprocess.run(cmd, shell=True, check=True)

def merge_files(dangerous_types, nondangerous_types):
    output_file = "merged.fasta"
    with open(output_file, 'w') as outfile:
        # First write the dangerous types
        for hpv_type in dangerous_types:
            filename = f"HPV_{hpv_type}_unique.fasta"
            with open(filename, 'r') as infile:
                outfile.write(infile.read())
        # Then write the non-dangerous types
        for hpv_type in nondangerous_types:
            filename = f"HPV_{hpv_type}_unique.fasta"
            with open(filename, 'r') as infile:
                outfile.write(infile.read())
    
    return output_file

def run_mafft(input_file, output_file):
    cmd = f"mafft --auto {input_file} > {output_file}" # align sequences with MAFFT
    subprocess.run(cmd, shell=True, check=True)

def calculate_mismatches_with_dangerous(probe, seq):
    sum = 0
    for a, b in zip(probe, seq):
        if a != b:
            sum += 1

    return sum

def calculate_mismatches_with_non_dangerous(probe, seq):
    sum = 0
    for a, b in zip(probe, seq):
        if a != b:
            sum += 1

    return sum

def select_probes(aligned_file, probe_length_range=range(25, 35), max_mismatches_dangerous=2, min_mismatches_nondangerous=3, probe_area_length=60):
    aligned_sequences = list(SeqIO.parse(aligned_file, 'fasta'))
    high_risk_seqs = [[] for i in range(36)]
    low_risk_seqs = [[] for i in range(45)]

    # Split sequences into high and low risk
    for i, record in enumerate(aligned_sequences):
        if i < 15:
            high_risk_seqs[0].append(str(record.seq).upper())
        elif i < 39:
            high_risk_seqs[1].append(str(record.seq).upper())
        elif i < 52:
            high_risk_seqs[2].append(str(record.seq).upper())
        elif i < 65:
            high_risk_seqs[3].append(str(record.seq).upper())
        elif i < 75:
            high_risk_seqs[4].append(str(record.seq).upper()) #all of the high risk sequences
        elif i < 102:
            low_risk_seqs[5].append(str(record.seq).upper())
        elif i < 118:
            low_risk_seqs[6].append(str(record.seq).upper())
        elif i < 123:
            low_risk_seqs[7].append(str(record.seq).upper())
        elif i < 129:
            low_risk_seqs[8].append(str(record.seq).upper())
        elif i < 137:
            low_risk_seqs[9].append(str(record.seq).upper())
        else:
            low_risk_seqs[10].append(str(record.seq).upper())

    selected_probes = set()
    found = False # flag to stop searching for probes once one is found in one type of virus

    for hr_seq_type in high_risk_seqs:
        found = False
        for hr_seq in hr_seq_type:
            for start in range(len(hr_seq) - max(probe_length_range) + 1):
                for length in probe_length_range:
                    probe_candidate = hr_seq[start:start+length]

                    if found:   # stop searching for probes once one is found in one type of virus
                        continue

                    if '-' in probe_candidate:  # skip probes with gaps
                        continue

                    if 'N' in probe_candidate:  # skip probes with unknown bases
                        continue

                    if len(probe_candidate) not in probe_length_range:  # skip probes with wrong length
                        continue

                    # Check if probe is similar to high risk sequences
                    for i, seq_type in enumerate(high_risk_seqs):
                        match = False   # flag to check if probe is similar to any sequence in the type
                        for seq in seq_type:
                            if calculate_mismatches_with_dangerous(probe_candidate, seq[start:start+length]) <= max_mismatches_dangerous:
                                match = True
                        if not match:
                            continue
                    
                    # Check if probe is unique to low risk sequences
                    for seq_type in low_risk_seqs:
                        for seq in seq_type:
                            if (hr_seq != seq):
                                if (calculate_mismatches_with_non_dangerous(probe_candidate, seq[start:start+length]) < min_mismatches_nondangerous):
                                    continue

                    selected_probes.add(probe_candidate)
                    found = True

    return list(selected_probes)

email = "private@gmail.com"
dangerous_types = ["16", "18", "31", "33", "35"]
non_dangerous_types = ["6", "11", "40", "42", "43", "44"]

# Download sequences for each HPV type
for hpv_type in dangerous_types + non_dangerous_types:
    print(f"Downloading HPV type {hpv_type} sequences...")
    download_sequences(email, hpv_type)

# Run CD-HIT on each HPV type
for hpv_type in dangerous_types + non_dangerous_types:
    print(f"Running CD-HIT on HPV type {hpv_type}...")
    run_cd_hit_wsl(f"HPV_{hpv_type}.fasta", f"HPV_{hpv_type}_unique.fasta")

# Merge all files into one
merged_file = merge_files(dangerous_types, non_dangerous_types)

# Run MAFFT on the merged file
run_mafft(merged_file, "aligned.fasta")

# Select probes
selected_probes = select_probes("aligned.fasta")

print("Selected Probes:")

for i, probe in enumerate(selected_probes, start=1):
    print(f"Probe {i}: {probe}")
