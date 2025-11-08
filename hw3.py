import sys

#取序列資料
def parse_fasta(path):
    """
    Parses a FASTA file.
    Assumes the file contains exactly two sequences.
    Returns a list of (header, sequence) tuples.
    """
    sequences = []
    header = ""
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:  # If we have a sequence, save it
                    sequences.append((header, seq))
                header = line[1:]  # Remove '>'
                seq = ""
            else:
                seq += line
        if header and seq:  # Add the last sequence
            sequences.append((header, seq))
    return sequences

#讀分數矩陣
def parse_score_matrix(path):
    """
    Parses a substitution matrix file (like PAM or BLOSUM).
    Returns a nested dictionary format: matrix[acid1][acid2] -> score
    """
    matrix = {}
    amino_acids = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue  # Skip comments
            
            parts = line.split()
            if not parts:
                continue

            if not amino_acids:
                # This is the header row
                amino_acids = parts
            else:
                # This is a score row
                row_acid = parts[0]
                matrix[row_acid] = {}
                scores = parts[1:]
                for i, col_acid in enumerate(amino_acids):
                    try:
                        matrix[row_acid][col_acid] = int(scores[i])
                    except (IndexError, ValueError):
                        print(f"Error parsing score file line: {line}", file=sys.stderr)
                        sys.exit(1)
    return matrix

#進行align_global - Needleman-Wunsch 全域比對
def align_global(seq1, seq2, score_matrix, gap):
    """
    Performs global alignment (Needleman-Wunsch).
    Returns one optimal alignment (aln1, aln2).
    """
    n = len(seq1)
    m = len(seq2)
    
    # Initialize DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[(0, 0) for _ in range(m + 1)] for _ in range(n + 1)] # (i-1, j-1), (i-1, j), (i, j-1)

    # Initialize first row and column
    for i in range(1, n + 1):
        dp[i][0] = dp[i-1][0] + gap
        trace[i][0] = (i - 1, 0) # From up
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j-1] + gap
        trace[0][j] = (0, j - 1) # From left

    # Fill the DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s1 = seq1[i-1]
            s2 = seq2[j-1]
            match_score = score_matrix[s1][s2]
            
            diag = dp[i-1][j-1] + match_score
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            
            scores = [diag, up, left]
            best_score = max(scores)
            dp[i][j] = best_score
            
            # Store traceback (prioritize diag, then up, then left)
            if best_score == diag:
                trace[i][j] = (i - 1, j - 1) # Diag
            elif best_score == up:
                trace[i][j] = (i - 1, j) # Up
            else:
                trace[i][j] = (i, j - 1) # Left

    # Traceback
    aln1, aln2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        prev_i, prev_j = trace[i][j]
        
        if prev_i == i - 1 and prev_j == j - 1: # Diagonal
            aln1 = seq1[i-1] + aln1
            aln2 = seq2[j-1] + aln2
            i, j = prev_i, prev_j
        elif prev_i == i - 1 and prev_j == j: # Up
            aln1 = seq1[i-1] + aln1
            aln2 = '-' + aln2
            i, j = prev_i, prev_j
        elif prev_i == i and prev_j == j - 1: # Left
            aln1 = '-' + aln1
            aln2 = seq2[j-1] + aln2
            i, j = prev_i, prev_j
        else:
            # Should only happen at (0,0) if one string is empty
            if i > 0:
                aln1 = seq1[i-1] + aln1
                aln2 = '-' + aln2
                i -= 1
            elif j > 0:
                aln1 = '-' + aln1
                aln2 = seq2[j-1] + aln2
                j -= 1
            else:
                break
                
    return [(aln1, aln2)] # Return as list for consistency

#align_local 進行 Smith-Waterman 區域比對
def align_local(seq1, seq2, score_matrix, gap):
    """
    Performs local alignment (Smith-Waterman).
    Handles complex tie-breaking.
    Returns a sorted list of optimal (aln1, aln2) tuples.
    """
    n = len(seq1)
    m = len(seq2)
    
    # DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    # Traceback matrix: stores a list of source directions
    trace = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    
    max_score = 0
    max_score_cells = [] # List of (i, j) tuples

    # Fill the DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s1 = seq1[i-1]
            s2 = seq2[j-1]
            match_score = score_matrix[s1][s2]
            
            diag = dp[i-1][j-1] + match_score
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            
            scores = [0, diag, up, left]
            
            best_score = max(scores)
            dp[i][j] = best_score
            
            # Store all directions that lead to the best score
            if best_score < 0:
                trace[i][j].append('start')
                #print('start : ',s1,', ',s2,' : '+'[',diag,', ', up,', ', left,']')
            else: # Only trace back if score is positive
                if best_score == diag:
                    trace[i][j].append('diag')
                    #print('diag  : ',s1,', ',s2,' : '+'[',diag,', ', up,', ', left,']')
                if best_score == up:
                    trace[i][j].append('up')
                    #print('up    : ',s1,', ',s2,' : '+'[',diag,', ', up,', ', left,']')
                if best_score == left:
                    trace[i][j].append('left')
                    #print('left  : ',s1,', ',s2,' : '+'[',diag,', ', up,', ', left,']')

            # Update max score and start points
            if best_score > max_score:
                max_score = best_score
                max_score_cells = [(i, j)]
            elif best_score == max_score and best_score > 0:
                max_score_cells.append((i, j))

    if not max_score_cells: # No positive alignment found
        return []

    # --- Multi-path Traceback ---
    # We use a set to store (aln1, aln2) tuples to avoid duplicates
    all_alignments = set()
    
    for start_i, start_j in max_score_cells:
        # Use a stack for iterative DFS to avoid recursion limits
        # State: (i, j, current_aln1, current_aln2)
        stack = [(start_i, start_j, "", "")]

        
        while stack:
            i, j, aln1, aln2 = stack.pop()
            
            '''
            # Base case: Hit a 'start' (score 0)
            if 'start' in trace[i][j] or dp[i][j] == 0:
                all_alignments.add((aln1, aln2))
                continue
            '''
            '''
            if dp[i][j] == 0 or 'start' in trace[i][j]:
                # include current cell before stopping
                all_alignments.add((seq1[i-1] + aln1, seq2[j-1] + aln2))
                continue
            ''' 
            '''
            if dp[i][j] == 0 or 'start' in trace[i][j]:
                all_alignments.add((aln1, aln2))
                continue
            ''' 
            
            
            for direction in trace[i][j]:
                
                if direction == 'diag':
                    stack.append((i - 1, j - 1, seq1[i-1] + aln1, seq2[j-1] + aln2))
                    print(direction,', ',seq1[i-1] + aln1,', ',seq2[j-1] + aln2)
                elif direction == 'up':
                    stack.append((i - 1, j, seq1[i-1] + aln1, '-' + aln2))
                elif direction == 'left':
                    stack.append((i, j - 1, '-' + aln1, seq2[j-1] + aln2))
            print(dp[i][j],', ',i,', ',j)
            if i == 0 or j == 0 or (dp[i][j] > dp[i-1][j-1]):
                all_alignments.add((aln1, aln2))
                continue
                
    if not all_alignments:
        return []
        
    # --- Filtering and Sorting ---
    
    # 1. Filter by max length
    # Note: len(aln1) == len(aln2) for all alignments
    max_len = 0
    for aln1, aln2 in all_alignments:
        if len(aln1) > max_len:
            max_len = len(aln1)
            
    max_len_alignments = [(a1, a2) for a1, a2 in all_alignments if len(a1) == max_len]
    
    # 2. Sort by aln1, then aln2
    # The set already handled duplicates
    sorted_alignments = sorted(max_len_alignments, key=lambda x: (x[0], x[1]))
    
    return sorted_alignments


def alignment(input_path, score_path, output_path, aln, gap):
    """
    Main alignment function.
    Reads inputs, calls the appropriate alignment algorithm,
    and writes the output FASTA file.
    """
    try:
        gap = int(gap)
    except ValueError:
        print(f"Error: Gap penalty '{gap}' must be an integer.", file=sys.stderr)
        sys.exit(1)

    try:
        seq_data = parse_fasta(input_path)
        print(seq_data)
        if len(seq_data) != 2:
            print(f"Error: Input FASTA '{input_path}' must contain exactly two sequences.", file=sys.stderr)
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing FASTA file '{input_path}': {e}", file=sys.stderr)
        sys.exit(1)

    try:
        score_matrix = parse_score_matrix(score_path)
    except FileNotFoundError:
        print(f"Error: Score file '{score_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing score file '{score_path}': {e}", file=sys.stderr)
        sys.exit(1)

    (header1, seq1) = seq_data[0]
    (header2, seq2) = seq_data[1]
    
    results = [] # List of (aln1, aln2) tuples

    try:
        if aln == "global":
            results = align_global(seq1, seq2, score_matrix, gap)
        elif aln == "local":
            results = align_local(seq1, seq2, score_matrix, gap)
        else:
            print(f"Error: Alignment type '{aln}' must be 'global' or 'local'.", file=sys.stderr)
            sys.exit(1)
        print(results)
    except KeyError as e:
        print(f"Error: Amino acid '{e.args[0]}' from input sequence not found in score matrix '{score_path}'.", file=sys.stderr)
        print("Please check that your sequences and score matrix are compatible.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during alignment: {e}", file=sys.stderr)
        sys.exit(1)

    # Write output
    try:
        with open(output_path, 'w') as f_out:
            if not results:
                # Write empty alignment if no result (e.g., local with no positive score)
                f_out.write(f">{header1}\n\n")
                f_out.write(f">{header2}\n\n")
            else:
                for aln1, aln2 in results:
                    f_out.write(f">{header1}\n{aln1}\n")
                    f_out.write(f">{header2}\n{aln2}\n")
    except IOError as e:
        print(f"Error writing to output file '{output_path}': {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    """
    Command line argument parser.
    Example usage:
    python alignment.py examples/test_global.fasta examples/pam250.txt examples/result_global.fasta global -10
    python alignment.py examples/test_local.fasta examples/pam100.txt examples/result_local.fasta local -10
    """
    if len(sys.argv) != 6:
        print("Usage: python alignment.py <input.fasta> <score.txt> <output.fasta> <global|local> <gap_score>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    score_file = sys.argv[2]
    output_file = sys.argv[3]
    aln_type = sys.argv[4]
    gap_penalty = sys.argv[5]
    
    alignment(input_file, score_file, output_file, aln_type, gap_penalty)
    
    # print(f"Alignment complete. Results written to {output_file}")