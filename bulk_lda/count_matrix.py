from operator import delitem
import pysam as ps
import pybedtools
import pandas as pd
import os
import os.path
import numpy as np
from pdb import set_trace
import pyBigWig

from .constants import MM_HEADER, MTX_SUFFIX, CELLS_SUFFIX, REGIONS_SUFFIX

def make_count_matrix(data: dict, output: str, merged_bed: str = "../data/merged_bed.bed", type = "dummy", norm = "rpkm") -> None:
    """Makes a count matrix given pairs of suitable peak calls and read data
    along with a normalisation method.

    Type can either be ["dummy", "bigwig", "bam"]


    Args:
        data (dict): Keys are peak files, values are their alignment files.
        output (str): Prefix to append .mtx .cells .regions
    """

    # First thing to do is to merge all of the peaks together and record
    # which peaks correspond to what. 
    # Read in all the keys as pybedtools instances
    corrected_data = {}
    for peak, bam in data.items():
        try:
            new_peak = append_file_names_to_bed(peak)
            corrected_data[new_peak] = bam
        except pd.errors.EmptyDataError:
            # It's okay if there's an empty peak file
            pass
 
    if not os.path.isfile(merged_bed):
        # Shove all the bed files together into a single one 
        # so that we can merge it 
        with open(merged_bed, 'w') as o:
            for fname in corrected_data.keys():
                with open(fname) as ifile:
                    for line in ifile:
                        o.write(line)
    
    merge = pybedtools.BedTool(merged_bed).sort().merge(c=11,o = "collapse").to_dataframe()

    if type.lower() == "dummy": 
        cell_reference = {k: i for i, k in enumerate(corrected_data.keys())}
        entries = []
        for peak_idx, row in merge.iterrows():
            cells = row['name'].split(",")
            for cell in cells:
                cell_idx = cell_reference[cell]
                reads = 1 
                entries.append([peak_idx+1, cell_idx+1, reads])
        
        final = np.array(entries)
    elif type.lower == "bigwig":
        tracks = {k: pyBigWig.open(v) for k, v in corrected_data.items()} 

        # Make sure they are all valid coverage tracks
        assert all([v.isBigWig() for k, v in tracks])

        cell_reference = {k: i for i, k in enumerate(corrected_data.keys())}
        entries = []
        for peak_idx, row in merge.iterrows():
            cells = row['name'].split(",")
            for cell in cells:
                cell_idx = cell_reference[cell]
                reads = len([i for i in tracks[cell].stats(row["chrom"], int(row["start"]), int(row["end"]))])
                peak_length = int(row["end"]) - int(row["start"])
                entries.append([peak_idx+1, cell_idx+1, reads, peak_length / 1000])

        final = np.array(entries)



    else: 
        # Loop through each of the peak entries and 
        #   1. figure out which files we have to look at
        #   2. find the read coverage under the peak 
        #   3. Add the entry as (peak index, cell index, read coverage)
        #       TODO: Need to discretise the read coverage somehow.
        #       Just use a list as apparently it is much faster

        # Open the bam files
        bams = {k: ps.AlignmentFile(v, "rb") for k, v in corrected_data.items()}

        # Starting with 0 here
        cell_reference = {k: i for i, k in enumerate(corrected_data.keys())}
        entries = []
        for peak_idx, row in merge.iterrows():
            cells = row['name'].split(",")
            for cell in cells:
                cell_idx = cell_reference[cell]
                reads = len([i for i in bams[cell].fetch(row["chrom"], int(row["start"]), int(row["end"]))])
                peak_length = int(row["end"]) - int(row["start"])
                entries.append([peak_idx+1, cell_idx+1, reads, peak_length / 1000])


        # Don't keep the connections open longer than I need them
        for key, bam in bams.items():
            bam.close()

        entries_np = np.array(entries)

        # Find all of the individual entries and split them up
        first = True
        for key, value in cell_reference.items():
            # Find where the condition is met
            # Need to normalise by
            #       total reads in the library
            #       length of the peak
            #
            # but they need to be an integer at the end of the day, 
            # so there needs to be some kind of transformation back to this.
            # what about quintiles? or something
            
            
            subset = entries_np[entries_np[:, 1] == value+1]
            
            if norm.lower() == "rpkm":
                total_m_reads = np.sum(subset[:, 2]) / 1e6  
                subset[:, 2] = np.round(subset[:, 2] / subset[:, 3] / total_m_reads)
                subset = subset[:, 0:3]
            elif norm.lower() == "none":
                subset = subset[:, 0:3]


            if first: 
                final = subset 
                first = False
            else:
                final = np.vstack((final, subset))

    # Sort the resulting matrix by peak index
    final = final[final[:, 0].argsort()]

    if os.path.isfile(f"{output}{MTX_SUFFIX}"):
        os.remove(f"{output}{MTX_SUFFIX}")

    # Write out the file
    with open(f"{output}{MTX_SUFFIX}", 'a') as mtx:
        mtx.write(MM_HEADER + "\n")
        mtx.write(f"{merge.shape[0]} {len(data.keys())} {int(np.sum(final[:, 2]))}\n")
        np.savetxt(mtx, final.astype(int), fmt = '%i', delimiter = " ")

    with open(f"{output}{REGIONS_SUFFIX}", 'w') as peaks_out:
        for i, row in merge.iterrows():
            peaks_out.write(f"{row.chrom}:{row.start}-{row.end}\n")
    
    with open(f"{output}{CELLS_SUFFIX}", 'w') as cells_out:
        for cell in cell_reference.keys():
            cells_out.write(cell + "\n")



def append_file_names_to_bed(bed: str) -> str:
    """Add file name annotation to bed file so that it can be collapsed.

    Does this as a side effect and returns the corrected file name
    
    Args:
        bed (str): Path to bed file
    Returns:
        str: New file name
    
    Side effects:
        Creates a new bed file with an additional column
    """
    new_file_name = bed.replace('.narrowPeak', '.new.bed')
    if not os.path.isfile(new_file_name):
        b = pd.read_csv(bed, sep = "\t", header = None)
        b['file'] = new_file_name 
        b.to_csv(new_file_name, sep = "\t", header = False, index = False)
    return new_file_name
