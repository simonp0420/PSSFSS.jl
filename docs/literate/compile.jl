using Literate

function notebook_filter(str)
  re1 = r"(?<!`)``(?!`)"  # Two backquotes not preceded by nor followed by another
  re2 = r"\[\^(\d)\]:?" # e.g. "[^1]" or "[^1]:" while capturing the digit
  str |> x -> replace(x, re1 => "\$") |> x -> replace(x, re2 => s"$^\1$")
end

#=
function md_filter(str)
    re = r"(^using .+Plots)"  # using at start of line followed by Plots later on line
    str |> x -> replace(x, re => s"""#md ENV["GKSwstype"] = "100" # hide\n\1""")
end
=#

examples_list = ["symmetric_strip.jl", "resistive_square_patch.jl", 
                 "cross_on_dielectric_substrate.jl", 
                 "square_loop_absorber.jl", "band_pass_filter.jl",
                 "cpss_optimization.jl"]

# Adds examples as subsections to the Examples page:
function replace_unknowns(str)
    findstr = """EditURL = "<unknown>"""
    replace(str, findstr => """EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate""")
end

flist = ["manual.jl"]
for file in flist
    Literate.markdown(file, "../src", codefence=("```@repl manual" => "```"), credit=false,
    postprocess=replace_unknowns)
    Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
end

for file in examples_list
  Literate.markdown(file, ".", credit=false)
  Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
end

function postinclude(str)
  str = replace_unknowns(str)
  for file in examples_list
    mdfile = splitext(file)[1] * ".md"
    str *= replace_unknowns(read(mdfile, String))
    rm(mdfile)
  end
  str
end

Literate.markdown("examples.jl", "../src", postprocess=postinclude, credit=false)

