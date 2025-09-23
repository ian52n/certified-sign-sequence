using CSV, DataFrames, Statistics
using StatsPlots  # groupedbar comes from StatsPlots

function analyze_repeat_probabilities(csv_file::String, min_count::Int=0)
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
        
        for row in eachrow(filtered_data)
            pattern = String(row["A"])
            event = String(row["Event"])
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
    
    println("\nFound repeat patterns:")
    println("L patterns: $(length(l_probs)) patterns")
    println("R patterns: $(length(r_probs)) patterns")
    
    # Sort by pattern length
    l_sorted = sortperm(l_lengths)
    r_sorted = sortperm(r_lengths)
    
    println("\nL repeat probabilities:")
    for i in l_sorted
        println("P(L|$(l_names[i])) = $(round(l_probs[i], digits=4)) [count: $(l_counts[i])]")
    end
    
    println("\nR repeat probabilities:")
    for i in r_sorted
        println("P(R|$(r_names[i])) = $(round(r_probs[i], digits=4)) [count: $(r_counts[i])]")
    end
    
    # FIRST PLOT
    p1 = plot(title="Repeat probabilities vs pattern length",
              xlabel="Pattern length",
              ylabel="P(repeat | pattern)",
              legend=:bottomright,
              size=(800, 600))
    
    if length(l_probs) > 0
        scatter!(p1, l_lengths[l_sorted], l_probs[l_sorted], 
                label="L patterns", color=:blue, markersize=6, alpha=0.8)
        plot!(p1, l_lengths[l_sorted], l_probs[l_sorted], 
              color=:blue, linewidth=2, alpha=0.6)
    end
    
    if length(r_probs) > 0
        scatter!(p1, r_lengths[r_sorted], r_probs[r_sorted], 
                label="R patterns", color=:red, markersize=6, alpha=0.8)
        plot!(p1, r_lengths[r_sorted], r_probs[r_sorted], 
              color=:red, linewidth=2, alpha=0.6)
    end
    
    if length(l_probs) > 0
        for i in l_sorted
            annotate!(p1, l_lengths[i], l_probs[i] + 0.01, text(l_names[i], 8, :blue))
        end
    end
    
    if length(r_probs) > 0
        for i in r_sorted
            annotate!(p1, r_lengths[i], r_probs[i] - 0.01, text(r_names[i], 8, :red))
        end
    end
    
    # SECOND PLOT: merged LX/RX style with true log₂ scaling, min y-value = 8
    lengths_sorted = sort(unique([l_lengths; r_lengths]))
    l_counts_by_len = [begin
        idx = findfirst(x -> x == len, l_lengths)
        idx === nothing ? 0 : l_counts[idx]
    end for len in lengths_sorted]
    r_counts_by_len = [begin
        idx = findfirst(x -> x == len, r_lengths)
        idx === nothing ? 0 : r_counts[idx]
    end for len in lengths_sorted]

    combined_labels = ["$(repeat('L', len))X / $(repeat('R', len))X" for len in lengths_sorted]

    # Transform counts to log₂ values for plotting
    log2_l_counts = [count > 0 ? log(count) / log(2) : 0 for count in l_counts_by_len]
    log2_r_counts = [count > 0 ? log(count) / log(2) : 0 for count in r_counts_by_len]
    counts_matrix_log2 = hcat(log2_l_counts, log2_r_counts)

    # Create yticks: powers of 2 starting at 8
    min_exp = ceil(Int, log(8) / log(2))  # 3
    max_count = maximum(vcat(l_counts_by_len, r_counts_by_len))
    max_exp = ceil(Int, log(max_count) / log(2))
    yticks_vals = collect(min_exp:max_exp)  # log₂ values
    yticks_labels = [string(round(2.0^val)) for val in yticks_vals]

    p2 = groupedbar(
        combined_labels,
        counts_matrix_log2,
        bar_position=:dodge,
        title="Counts for event A in P(B | A) by pattern length (log₂ scale)",
        xlabel="Pattern (L vs R)",
        ylabel="Count(A) [log₂ scale]",
        size=(800, 600),
        color=[:blue :red],
        xticks=(1:length(combined_labels), combined_labels),
        xrotation=45,
        yticks=(yticks_vals, yticks_labels),
        ylim=(min_exp, max_exp)  # set minimum y-value to log₂(8)
    )
    
    combined_plot = plot(p1, p2, layout=(2,1), size=(900, 1100), margin=10Plots.mm)
    
    output_file = replace(csv_file, ".csv" => "_repeat_analysis.png")
    savefig(combined_plot, output_file)
    println("\nRepeat analysis plot saved to: $output_file")
    
    println("\nStatistical Summary:")
    println("==================")
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

if length(ARGS) >= 1
    csv_file = ARGS[1]
    min_count = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 0
    try
        l_probs, l_lengths, r_probs, r_lengths, matched_l, matched_r = analyze_repeat_probabilities(csv_file, min_count)
    catch e
        println("Error: $e")
        rethrow(e)
    end
else
    println("Usage: julia RepeatProbabilityAnalysis.jl analysis_results.csv [min_count]")
    println("Default min_count is 0")
    println("This will analyze how repeat probabilities change with pattern length")
end
