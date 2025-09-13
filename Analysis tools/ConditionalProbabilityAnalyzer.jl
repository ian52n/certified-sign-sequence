using CSV, DataFrames

function analyze_sequences(input_file::String, output_file::String)
    # Read sequences from file
    sequences = String[]
    
    if !isfile(input_file)
        error("File not found: $input_file")
    end
    
    open(input_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line)
                push!(sequences, line)
            end
        end
    end
    
    if isempty(sequences)
        error("No sequences found in file")
    end
    
    # Find maximum sequence length
    max_length = maximum(length.(sequences))
    println("Maximum sequence length: $max_length")
    
    # We'll compute statistics for n-grams up to a reasonable maximum
    # For very long sequences, limit to prevent exponential explosion
    max_ngram = min(max_length - 1, 10)  # Cap at 10-grams for performance
    
    if max_ngram < 1
        error("All sequences are too short (need at least length 2)")
    end
    
    println("Computing statistics for 1-gram through $(max_ngram)-gram")
    
    # Generate all possible n-gram patterns for each n
    function generate_patterns(n::Int)
        if n == 1
            return ["L", "R"]
        end
        patterns = String[]
        for i in 0:(2^n - 1)
            pattern = ""
            temp = i
            for j in 1:n
                if temp % 2 == 0
                    pattern = "L" * pattern
                else
                    pattern = "R" * pattern
                end
                temp ÷= 2
            end
            push!(patterns, pattern)
        end
        return patterns
    end
    
    # Initialize data structures for all n-grams
    # pattern_counts[n][pattern] = count of pattern that can be followed by another character
    pattern_counts = Dict{Int, Dict{String, Int}}()
    # transition_counts[n][context][next_char] = count of next_char following context
    transition_counts = Dict{Int, Dict{String, Dict{Char, Int}}}()
    
    for n in 1:max_ngram
        pattern_counts[n] = Dict{String, Int}()
        transition_counts[n] = Dict{String, Dict{Char, Int}}()
        
        # Initialize all possible patterns
        for pattern in generate_patterns(n)
            pattern_counts[n][pattern] = 0
            transition_counts[n][pattern] = Dict('L' => 0, 'R' => 0)
        end
    end
    
    # Process each sequence
    for seq in sequences
        seq_length = length(seq)
        
        # Process all n-grams for this sequence
        for n in 1:min(max_ngram, seq_length - 1)
            
            # CRITICAL: Count n-gram patterns only if they can be followed by another character
            # This means we count positions 1 through (seq_length - n) inclusive
            # For a sequence of length L, an n-gram at position i uses characters i through i+n-1
            # For this n-gram to be "followable", position i+n must exist (i+n ≤ L)
            # So we need i ≤ L - n, which gives us positions 1 through (L - n)
            
            for i in 1:(seq_length - n)
                pattern = seq[i:i+n-1]
                
                # Verify this pattern can be followed by checking if position i+n exists
                if i + n <= seq_length
                    if haskey(pattern_counts[n], pattern)
                        pattern_counts[n][pattern] += 1
                    end
                    
                    # Count the transition
                    next_char = seq[i+n]
                    if haskey(transition_counts[n], pattern)
                        transition_counts[n][pattern][next_char] += 1
                    end
                end
            end
        end
    end
    
    # Build results in the new format
    result_data = []
    
    # Process results for each n-gram level
    for n in 1:max_ngram
        println("Processing $(n)-grams...")
        
        patterns = generate_patterns(n)
        
        for pattern in sort(patterns)  # Sort for consistent output
            count_pattern = pattern_counts[n][pattern]
            
            # For each possible next character, add a row
            for next_char in ['L', 'R']
                total_transitions = sum(values(transition_counts[n][pattern]))
                
                if total_transitions > 0
                    transition_count = transition_counts[n][pattern][next_char]
                    prob = transition_count / total_transitions
                else
                    prob = 0.0
                end
                
                event_notation = string(next_char) * "|" * pattern
                
                push!(result_data, (
                    A = pattern,
                    CountA = count_pattern,
                    Event = event_notation,
                    ProbBA = prob
                ))
            end
        end
    end
    
    # Create DataFrame from the collected data
    A_col = [row.A for row in result_data]
    CountA_col = [row.CountA for row in result_data]
    Event_col = [row.Event for row in result_data]
    ProbBA_col = [row.ProbBA for row in result_data]
    
    results = DataFrame()
    results[!, "A"] = A_col
    results[!, "Count(A)"] = CountA_col
    results[!, "Event"] = Event_col
    results[!, "P(B|A)"] = ProbBA_col
    
    # Write to CSV
    CSV.write(output_file, results)
    
    println("Analysis complete. Results saved to: $output_file")
    println("Processed $(length(sequences)) sequences")
    println("Computed statistics for 1-gram through $(max_ngram)-gram")
    
    # Print summary statistics
    total_rows = length(result_data)
    unique_patterns = length(unique(A_col))
    println("Total rows in output: $total_rows")
    println("Unique patterns analyzed: $unique_patterns")
    
    # Show a few example rows
    if total_rows > 0
        println("\nFirst few rows of results:")
        println("A\t\tCount(A)\tEvent\t\tP(B|A)")
        for i in 1:min(6, total_rows)
            println("$(A_col[i])\t\t$(CountA_col[i])\t\t$(Event_col[i])\t\t$(round(ProbBA_col[i], digits=4))")
        end
    end
end

# Example usage:
if length(ARGS) >= 2
    input_file = ARGS[1]
    output_file = ARGS[2]
    try
        analyze_sequences(input_file, output_file)
    catch e
        println("Error: $e")
        rethrow(e)
    end
elseif length(ARGS) == 1
    input_file = ARGS[1]
    output_file = replace(input_file, ".txt" => "_analysis.csv")
    try
        analyze_sequences(input_file, output_file)
    catch e
        println("Error: $e")
        rethrow(e)
    end
else
    println("Usage: julia script.jl input_file.txt [output_file.csv]")
    println("If output file is not specified, it will be named input_file_analysis.csv")
end