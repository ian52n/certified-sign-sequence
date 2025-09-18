"""
    invert_word(word::String) -> String

Inverts a word made of 'L' and 'R' characters.
'L' becomes 'R' and 'R' becomes 'L'.
"""
function invert_word(word::String)
    return join(c == 'L' ? 'R' : 'L' for c in word)
end

"""
    verify_pattern(filename::String)

Reads a file line by line and checks if words come in pairs
where the second word is an inversion of the first.
"""
function verify_pattern(filename::String)
    # Check if the file exists before trying to read it
    if !isfile(filename)
        println("❌ Error: File '$filename' not found.")
        return
    end

    lines = readlines(filename)
    word_count = length(lines)
    
    # The pattern requires an even number of words to form pairs
    if isodd(word_count)
        println("❌ Pattern is invalid: The file contains $word_count words, which is an odd number.")
        return
    end

    println("🔎 Checking $word_count words in '$filename' for the inversion pattern...")

    # Loop through the lines by steps of 2 to process pairs
    for i in 1:2:word_count
        word1 = lines[i]
        word2 = lines[i+1] # The next word in the pair

        # Generate the expected inverted version of the first word
        expected_word2 = invert_word(word1)

        # Check if the actual second word matches the expected inversion
        if word2 != expected_word2
            println("\n❌ Pattern broken at lines $i and $(i+1):")
            println("   Word 1:           '$word1'")
            println("   Expected Word 2:  '$expected_word2'")
            println("   Actual Word 2:    '$word2'")
            return # Exit function as soon as a mismatch is found
        end
    end

    # If the loop completes without finding any mismatches
    println("\n✅ Success! The L/R inversion pattern holds true for all words.")
end

# --- Main execution ---
# Check if the ARGS array is empty. ARGS holds command-line arguments.
if isempty(ARGS)
    println("Usage: julia $(@__FILE__) <filename>")
else
    # Use the first command-line argument as the filename
    filename = ARGS[1]
    verify_pattern(filename)
end