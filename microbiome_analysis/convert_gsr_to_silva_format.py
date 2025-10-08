#!/usr/bin/env python3
"""
Script to convert GSR-DB fasta file to SILVA format by combining sequence IDs with taxonomic information.

Input files:
- GSR-DB_V3-V4_cluster-1_seqs.fasta: FASTA file with simple sequence IDs
- GSR-DB_V3-V4_cluster-1_taxa.txt: Tab-delimited file with Feature ID and taxonomic lineage

Output:
- GSR-DB_V3-V4_cluster-1_silva_format.fasta.gz: Compressed FASTA file with SILVA-style headers

SILVA format: >SequenceID TaxonomicLineage;OrganismName
"""

import gzip
import sys
import argparse
from pathlib import Path


def parse_taxonomy_file(taxa_file):
    """
    Parse the taxonomy file and return a dictionary mapping sequence IDs to taxonomic lineages.
    
    Args:
        taxa_file (str): Path to the tab-delimited taxonomy file
        
    Returns:
        dict: Dictionary mapping sequence IDs to taxonomic lineages
    """
    taxonomy_dict = {}
    
    print(f"Reading taxonomy file: {taxa_file}")
    
    with open(taxa_file, 'r') as f:
        # Skip header line
        next(f)
        
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"Warning: Line {line_num} has unexpected format: {line}")
                continue
                
            seq_id = parts[0].strip()
            taxon = parts[1].strip()
            
            # Convert GSR-DB format to SILVA format
            # GSR-DB: k__Bacteria; p__Actinobacteria; c__Actinomycetia; ...
            # SILVA: Bacteria;Actinobacteria;Actinomycetia; ...
            silva_taxon = convert_taxon_format(taxon)
            
            taxonomy_dict[seq_id] = silva_taxon
    
    print(f"Loaded {len(taxonomy_dict)} taxonomic entries")
    return taxonomy_dict


def convert_taxon_format(taxon):
    """
    Convert GSR-DB taxonomic format to SILVA format.
    Returns only 7 taxonomic levels: Kingdom through Species.
    
    Args:
        taxon (str): Taxonomic string in GSR-DB format
        
    Returns:
        str: Taxonomic string in SILVA format (7 levels only)
    """
    # Split by semicolon and remove prefixes (k__, p__, c__, etc.)
    parts = taxon.split(';')
    silva_parts = []
    
    for part in parts:
        part = part.strip()
        if '__' in part:
            # Remove prefix (k__, p__, c__, o__, f__, g__, s__)
            silva_part = part.split('__', 1)[1]
            # Handle multiple options separated by hyphens
            if '-' in silva_part and 'Unknown' not in silva_part:
                # Take the first option when multiple are available
                silva_part = silva_part.split('-')[0]
            silva_parts.append(silva_part)
        else:
            silva_parts.append(part)
    
    # Limit to 7 taxonomic levels: Kingdom, Phylum, Class, Order, Family, Genus, Species
    # SILVA format should only have these 7 levels, not 8
    if len(silva_parts) > 7:
        silva_parts = silva_parts[:7]
    
    # Join with semicolons (no additional organism name)
    silva_taxon = ';'.join(silva_parts)
    
    return silva_taxon


def convert_fasta_to_silva_format(fasta_file, taxonomy_dict, output_file):
    """
    Convert FASTA file to SILVA format by adding taxonomic information to headers.
    
    Args:
        fasta_file (str): Path to input FASTA file
        taxonomy_dict (dict): Dictionary mapping sequence IDs to taxonomic lineages
        output_file (str): Path to output compressed FASTA file
    """
    print(f"Converting FASTA file: {fasta_file}")
    print(f"Output file: {output_file}")
    
    sequences_processed = 0
    sequences_matched = 0
    sequences_missing_taxonomy = 0
    
    with open(fasta_file, 'r') as infile, gzip.open(output_file, 'wt') as outfile:
        current_seq_id = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # Process previous sequence if exists
                if current_seq_id and current_sequence:
                    sequences_processed += 1
                    
                    if current_seq_id in taxonomy_dict:
                        sequences_matched += 1
                        silva_header = f">{current_seq_id} {taxonomy_dict[current_seq_id]}"
                        outfile.write(silva_header + '\n')
                        outfile.write('\n'.join(current_sequence) + '\n')
                    else:
                        sequences_missing_taxonomy += 1
                        print(f"Warning: No taxonomy found for sequence {current_seq_id}")
                        # Still write the sequence with minimal taxonomic info
                        silva_header = f">{current_seq_id} Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
                        outfile.write(silva_header + '\n')
                        outfile.write('\n'.join(current_sequence) + '\n')
                
                # Start new sequence
                current_seq_id = line[1:]  # Remove '>' prefix
                current_sequence = []
            else:
                # Add sequence line
                current_sequence.append(line)
        
        # Process last sequence
        if current_seq_id and current_sequence:
            sequences_processed += 1
            
            if current_seq_id in taxonomy_dict:
                sequences_matched += 1
                silva_header = f">{current_seq_id} {taxonomy_dict[current_seq_id]}"
                outfile.write(silva_header + '\n')
                outfile.write('\n'.join(current_sequence) + '\n')
            else:
                sequences_missing_taxonomy += 1
                print(f"Warning: No taxonomy found for sequence {current_seq_id}")
                silva_header = f">{current_seq_id} Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
                outfile.write(silva_header + '\n')
                outfile.write('\n'.join(current_sequence) + '\n')
    
    print(f"\nConversion completed:")
    print(f"  Total sequences processed: {sequences_processed}")
    print(f"  Sequences with taxonomy: {sequences_matched}")
    print(f"  Sequences missing taxonomy: {sequences_missing_taxonomy}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GSR-DB FASTA file to SILVA format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_gsr_to_silva_format.py \\
    --fasta references/GSR-DB_V3-V4_cluster-1_seqs.fasta \\
    --taxa references/GSR-DB_V3-V4_cluster-1_taxa.txt \\
    --output references/GSR-DB_V3-V4_cluster-1_silva_format.fasta.gz
        """
    )
    
    parser.add_argument(
        '--fasta', 
        required=True,
        help='Input FASTA file with sequence IDs'
    )
    
    parser.add_argument(
        '--taxa',
        required=True, 
        help='Input taxonomy file (tab-delimited with Feature ID and Taxon columns)'
    )
    
    parser.add_argument(
        '--output',
        required=True,
        help='Output compressed FASTA file in SILVA format'
    )
    
    args = parser.parse_args()
    
    # Validate input files
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}")
        sys.exit(1)
        
    if not Path(args.taxa).exists():
        print(f"Error: Taxonomy file not found: {args.taxa}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Parse taxonomy file
        taxonomy_dict = parse_taxonomy_file(args.taxa)
        
        # Convert FASTA file
        convert_fasta_to_silva_format(args.fasta, taxonomy_dict, args.output)
        
        print(f"\nSuccess! Converted FASTA file saved as: {args.output}")
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
