using Statistics, HypothesisTests

# Count total alphabetic letters and R's in file
function count_letters_and_Rs(filename::String; case_insensitive::Bool=true)
    text = read(filename, String)
    letters = filter(isletter, text)               # keep only alphabetic chars
    n = length(letters)
    if case_insensitive
        countR = count(c -> lowercase(string(c)) == "r", letters)
    else
        countR = count(c -> c == 'R', letters)
    end
    return n, countR
end

# --- Main ---
if length(ARGS) < 1
    error("Usage: julia Bernoullitest.jl <filename> [--case-sensitive]")
end

filename = ARGS[1]
cs_flag = any(a -> a == "--case-sensitive", ARGS)  # optional flag to count only 'R'

n, countR = count_letters_and_Rs(filename; case_insensitive = !cs_flag)

if n == 0
    error("No alphabetic letters found in file: $filename")
end

countL = n - countR

# Build test and get two-sided p-value
test = BinomialTest(countR, n, 0.5)   # by default this is two-sided
pval = pvalue(test)

println("File: $filename")
println("Total letters (n) = $n")
println("R = $countR, L = $countL")
println("Binomial test (H0: p=0.5) two-sided p-value = $pval")
