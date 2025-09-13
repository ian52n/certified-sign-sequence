using CSV, DataFrames, Plots, Statistics

function check_lr_symmetry(csv_file::String)
    # Read the analysis results
    data = CSV.read(csv_file, DataFrame)
    
    println("Loaded $(nrow(data)) rows from $csv_file")
    
    # Count total occurrences of valid strings (sum of Count(A) where Count(A) >= 100)
    total_occurrences = sum(data[data[!, "Count(A)"] .>= 100, "Count(A)"])
    println("Total occurrences of strings with Count(A) >= 100: $total_occurrences")
    
    # Function to flip L's and R's in a string
    function flip_lr(s::String)
        return replace(replace(s, 'L' => 'X'), 'R' => 'L', 'X' => 'R')
    end
    
    # Create pairs for comparison
    original_probs = Float64[]
    flipped_probs = Float64[]
    pattern_labels = String[]
    count_weights = Int64[]  # Add weights for each comparison
    
    # For each row in the original data
    for i in 1:nrow(data)
        original_pattern = String(data[i, "A"])
        original_event = String(data[i, "Event"]) 
        original_prob = data[i, "P(B|A)"]
        original_count = data[i, "Count(A)"]
        
        # Create the flipped version
        flipped_pattern = flip_lr(original_pattern)
        flipped_event = flip_lr(original_event)
        
        # Find the matching flipped row in the data
        matching_rows = data[(data[!, "A"] .== flipped_pattern) .& (data[!, "Event"] .== flipped_event), :]
        
        if nrow(matching_rows) == 1
            flipped_count = matching_rows[1, "Count(A)"]
            flipped_prob = matching_rows[1, "P(B|A)"]
            
            # Only include if both original and flipped have sufficient counts
            if flipped_count >= 100 && original_count >= 100
                push!(original_probs, original_prob)
                push!(flipped_probs, flipped_prob)
                push!(pattern_labels, original_event)
                push!(count_weights, original_count)  # Use original count as weight
            end
        end
    end
    
    println("Found $(length(original_probs)) symmetric pairs")
    
    # Create scatter plot
    p1 = scatter(original_probs, flipped_probs, 
                 xlabel="Original Probability P(B|A)", 
                 ylabel="L/R Flipped Probability P(B'|A')",
                 title="L/R Symmetry Check",
                 alpha=0.7,
                 markersize=3,
                 legend=false)
    
    # Add perfect symmetry line (y=x)
    plot!(p1, [0, 1], [0, 1], 
          linestyle=:dash, 
          linewidth=2, 
          color=:red,
          label="Perfect Symmetry (y=x)")
    
    # Calculate and display correlation
    correlation = cor(original_probs, flipped_probs)
    annotate!(p1, 0.1, 0.9, text("Correlation: $(round(correlation, digits=4))", 10))
    
    # Calculate mean absolute difference (weighted by counts)
    differences = original_probs .- flipped_probs
    mad = sum(abs.(differences) .* count_weights) / sum(count_weights)
    annotate!(p1, 0.1, 0.85, text("Weighted Mean Abs Diff: $(round(mad, digits=4))", 10))
    
    # Create histogram of differences (weighted by counts)
    # Create weighted histogram by repeating each difference count_weights[i] times
    weighted_differences = Float64[]
    for i in 1:length(differences)
        append!(weighted_differences, fill(differences[i], count_weights[i]))
    end
    
    p2 = histogram(weighted_differences,
                   bins=20,
                   xlabel="Difference (Original - Flipped)",
                   ylabel="Weighted Count",
                   title="Weighted Distribution of Probability Differences",
                   alpha=0.7,
                   legend=false)
    
    # Add vertical line at zero
    vline!(p2, [0], linestyle=:dash, linewidth=2, color=:red)
    
    # Create combined plot with better margins
    combined_plot = plot(p1, p2, layout=(1,2), size=(1200, 500), 
                        margin=15Plots.mm, # Add margins around the entire plot
                        left_margin=10Plots.mm, right_margin=10Plots.mm,
                        top_margin=5Plots.mm, bottom_margin=10Plots.mm)
    
    # Save the plot
    output_file = replace(csv_file, ".csv" => "_symmetry_check.png")
    savefig(combined_plot, output_file)
    println("Symmetry check plot saved to: $output_file")
    
    # Print some summary statistics
    println("\nSymmetry Analysis Summary:")
    println("========================")
    println("Correlation coefficient: $(round(correlation, digits=4))")
    println("Weighted mean absolute difference: $(round(mad, digits=4))")
    weighted_std = sqrt(sum(((differences .- sum(differences .* count_weights) / sum(count_weights)).^2) .* count_weights) / sum(count_weights))
    println("Weighted standard deviation of differences: $(round(weighted_std, digits=4))")
    println("Maximum difference: $(round(maximum(abs.(differences)), digits=4))")
    println("Total weighted observations: $(sum(count_weights))")
    
    # Find and report the most asymmetric patterns (top 5 by absolute difference, showing weights)
    abs_diffs = abs.(differences)
    sorted_indices = sortperm(abs_diffs, rev=true)
    
    println("\nMost asymmetric patterns (top 5 by absolute difference):")
    for i in 1:min(5, length(sorted_indices))
        idx = sorted_indices[i]
        println("$(pattern_labels[idx]) [count: $(count_weights[idx])]: $(round(original_probs[idx], digits=4)) vs $(round(flipped_probs[idx], digits=4)) (diff: $(round(differences[idx], digits=4)))")
    end
    
    return correlation, mad, differences
end

# Example usage:
if length(ARGS) >= 1
    csv_file = ARGS[1]
    try
        correlation, mad, differences = check_lr_symmetry(csv_file)
    catch e
        println("Error: $e")
        rethrow(e)
    end
else
    println("Usage: julia symmetry_check.jl analysis_results.csv")
    println("This will create a visualization showing L/R symmetry in your probability data")
end