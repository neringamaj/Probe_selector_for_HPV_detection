<h1>HPV Probe Selection System for PCR</h1>
<h2>Overview</h2>
<p>
    This project develops a program for the selective detection of high-risk Human Papillomavirus (HPV) types that are known to cause cervical cancer, using Polymerase Chain Reaction (PCR) for DNA amplification. The system is tailored to identify HPV types 16, 18, 31, 33, and 35 while ensuring specificity by not signaling the presence of non-dangerous HPV types 6, 11, 40, 42, 43, and 44.
</p>

<h2>Methodology</h2>
<p>
    The program automates the selection of DNA probes within the L1 gene region, optimizing for the smallest set of probes that cover all high-risk HPV types. The selection criteria for the probes are as follows:
</p>
<ul>
    <li>Probe length: 30±5 base pairs (bp).</li>
    <li>Match specificity: At least one probe matches each high-risk HPV sequence with no more than 2 mismatches.</li>
    <li>Non-dangerous HPV specificity: At least 3 mismatches with all non-dangerous HPV sequences.</li>
</ul>
<p>
    The process involves:
</p>
<ol>
    <li>Downloading both dangerous and non-dangerous HPV sequences in FASTA format.</li>
    <li>Utilizing the CD-HIT program to eliminate identical sequences.</li>
    <li>Merging the sequences into a single FASTA file and conducting a comprehensive comparison using MAFFT.</li>
    <li>Selecting a suitable probe system based on the above criteria.</li>
</ol>

<h2>Requirements</h2>
<p>
    <ul>
        <li>CD-HIT</li>
        <li>MAFFT</li>
        <li>Python</li>
        <li>BioPython (for sequence handling and BLAST queries)</li>
    </ul>
</p>
