#!/usr/bin/env julia
# extract_signs.jl
# Usage: julia extract_signs.jl [input_dir] [output_file]
# Default input_dir = "." ; default output_file = "./combined_sign_sequences.txt"

# Extracts sign sequences from the bottom of .txt files and puts them
# all in a single .txt file

using Printf

const MARKER = "Sign sequences separated by newlines:"

function extract_from_file(path::AbstractString; marker::AbstractString=MARKER)
    text = read(path, String)
    lines = split(text, '\n')
    idx = findfirst(line -> strip(line) == marker, lines)
    if idx === nothing
        return nothing
    end
    if idx >= length(lines)
        return ""  # marker is last line, nothing after it
    end
    return join(lines[idx+1:end], "\n")
end

function main(argv)
    input_dir = length(argv) >= 1 ? argv[1] : "."
    output_path = length(argv) >= 2 ? argv[2] : joinpath(".", "combined_sign_sequences.txt")

    files = sort(filter(f -> lowercase(f) |> endswith(".txt") && isfile(joinpath(input_dir, f)), readdir(input_dir)))
    n_written = 0
    open(output_path, "w") do out
        first_written = false
        prev_ended_nl = false
        for f in files
            path = joinpath(input_dir, f)
            chunk = extract_from_file(path)
            if chunk === nothing
                @warn "marker not found; skipping file" file=path
                continue
            end
            if isempty(chunk)
                continue
            end
            if first_written && !prev_ended_nl
                write(out, "\n")
            end
            write(out, chunk)
            prev_ended_nl = endswith(chunk, "\n")
            first_written = true
            n_written += 1
        end
    end

    @printf("Extracted sign sequences from %d file(s) into %s\n", n_written, output_path)
end

main(ARGS)
