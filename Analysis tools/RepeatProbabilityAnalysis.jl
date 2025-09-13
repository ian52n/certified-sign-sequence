using CSV, DataFrames, Plots, Statistics

function analyze_repeat_probabilities(csv_file::String, min_count::Int=30)
    # Read the analysis results
    data = CSV.read(csv_file, DataFrame)
    
    println("Loaded $(nrow(data)) rows from $csv_file")
    println("Using minimum count threshold: $min_count")
    
    # Filter data for sufficient counts
    filtered_data = data[data[!, "Count(A)"] .>= min_count, :]
    println("After filtering: $(nrow(filtered_data)) rows")
    
    # Function to extract repeat probabilities
    function extract_repeat_probs(letter::Char)
        repeat_probs = Float64[]
        pattern_lengths = Int[]
        pattern_names = String[]
        counts = Int[]
        
        # Find all patterns that are repeats of the given letter
        for row in eachrow(filtered_data)
            pattern = String(row["A"])
            event = String(row["Event"])
            
            # Check if pattern is all the same letter and event is asking for that same letter
            if all(c -> c == letter, pattern) && event == string(letter, "|", pattern)
                push!(repeat_probs, row["P(B|A)"])
                push!(pattern_lengths, length(pattern))
                push!(pattern_names, pattern)
                push!(counts, row["Count(A)"])
            end
        end
        
        return repeat_probs, pattern_lengths, pattern_names, counts
    end
    
    # Extract data for L and R
    l_probs, l_lengths, l_names, l_counts = extract_repeat_probs('L')
    r_probs, r_lengths, r_names, r_counts = extract_repeat_probs('R')
    
    println("\\nFound repeat patterns:")
    println("L patterns: $(length(l_probs)) patterns")
    println("R patterns: $(length(r_probs)) patterns")
    
    # Sort by pattern length
    l_sorted = sortperm(l_lengths)
    r_sorted = sortperm(r_lengths)
    
    # Print the sequences
    println("\\nL repeat probabilities:")
    for i in l_sorted
        println("P(L|$(l_names[i])) = $(round(l_probs[i], digits=4)) [count: $(l_counts[i])]")
    end
    
    println("\\nR repeat probabilities:")
    for i in r_sorted
        println("P(R|$(r_names[i])) = $(round(r_probs[i], digits=4)) [count: $(r_counts[i])]")
    end
    
    # Create visualization
    p1 = plot(title="Repeat Probabilities vs Pattern Length",
              xlabel="Pattern Length",
              ylabel="P(repeat | pattern)",
              legend=:bottomright,
              size=(800, 600))
    
    # Plot L patterns
    if length(l_probs) > 0
        scatter!(p1, l_lengths[l_sorted], l_probs[l_sorted], 
                label="L patterns", 
                color=:blue, 
                markersize=6,
                alpha=0.8)
        # Connect points with lines
        plot!(p1, l_lengths[l_sorted], l_probs[l_sorted], 
              color=:blue, 
              linewidth=2,
              alpha=0.6)
    end
    
    # Plot R patterns  
    if length(r_probs) > 0
        scatter!(p1, r_lengths[r_sorted], r_probs[r_sorted], 
                label="R patterns", 
                color=:red, 
                markersize=6,
                alpha=0.8)
        # Connect points with lines
        plot!(p1, r_lengths[r_sorted], r_probs[r_sorted], 
              color=:red, 
              linewidth=2,
              alpha=0.6)
    end
    
    # Add pattern labels
    if length(l_probs) > 0
        for i in l_sorted
            annotate!(p1, l_lengths[i], l_probs[i] + 0.01, 
                     text(l_names[i], 8, :blue))
        end
    end
    
    if length(r_probs) > 0
        for i in r_sorted
            annotate!(p1, r_lengths[i], r_probs[i] - 0.01, 
                     text(r_names[i], 8, :red))
        end
    end
    
    # Create a comparison plot showing L vs R symmetry
    p2 = plot(title="L/R Symmetry in Repeat Probabilities",
              xlabel="L Repeat Probability",
              ylabel="R Repeat Probability",
              legend=:bottomright,
              size=(600, 600))
    
    # Match L and R patterns of same length
    matched_l = Float64[]
    matched_r = Float64[]
    matched_lengths = Int[]
    
    for len in unique([l_lengths; r_lengths])
        l_idx = findfirst(x -> x == len, l_lengths)
        r_idx = findfirst(x -> x == len, r_lengths)
        
        if l_idx !== nothing && r_idx !== nothing
            push!(matched_l, l_probs[l_idx])
            push!(matched_r, r_probs[r_idx])
            push!(matched_lengths, len)
        end
    end
    
    if length(matched_l) > 0
        scatter!(p2, matched_l, matched_r, 
                label="Matched lengths",
                color=:purple,
                markersize=8,
                alpha=0.8)
        
        # Add y=x line for perfect symmetry
        plot!(p2, [0, 1], [0, 1], 
              linestyle=:dash,
              color=:black,
              linewidth=2,
              label="Perfect symmetry")
        
        # Add length labels
        for i in 1:length(matched_l)
            annotate!(p2, matched_l[i], matched_r[i], 
                     text(string(matched_lengths[i]), 10, :purple))
        end
        
        # Calculate correlation
        if length(matched_l) > 1
            corr = cor(matched_l, matched_r)
            annotate!(p2, 0.1, 0.9, 
                     text("Correlation: $(round(corr, digits=4))", 10))
        end
    end
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 1000),
                        margin=10Plots.mm)
    
    # Save the plot
    output_file = replace(csv_file, ".csv" => "_repeat_analysis.png")
    savefig(combined_plot, output_file)
    println("\\nRepeat analysis plot saved to: $output_file")
    
    # Statistical summary
    println("\\nStatistical Summary:")
    println("==================")
    if length(matched_l) > 0
        println("Matched pattern pairs: $(length(matched_l))")
        if length(matched_l) > 1
            println("L/R correlation: $(round(cor(matched_l, matched_r), digits=4))")
            println("Mean L repeat prob: $(round(mean(matched_l), digits=4))")
            println("Mean R repeat prob: $(round(mean(matched_r), digits=4))")
            println("Mean absolute difference: $(round(mean(abs.(matched_l .- matched_r)), digits=4))")
        end
    else
        println("No matched L/R patterns found at this count threshold")
    end
    
    return l_probs, l_lengths, r_probs, r_lengths, matched_l, matched_r
end

# Example usage:
if length(ARGS) >= 1
    csv_file = ARGS[1]
    min_count = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 30
    
    try
        l_probs, l_lengths, r_probs, r_lengths, matched_l, matched_r = analyze_repeat_probabilities(csv_file, min_count)
    catch e
        println("Error: $e")
        rethrow(e)
    end
else
    println("Usage: julia repeat_analysis.jl analysis_results.csv [min_count]")
    println("Default min_count is 30")
    println("This will analyze how repeat probabilities change with pattern length")
end