#using Literate

function notebook_filter(str)
  re1 = r"(?<!`)``(?!`)"  # Two backquotes not preceded by nor followed by another
  re2 = r"\[\^(\d)\]:?" # e.g. "[^1]" or "[^1]:" while capturing the digit
  str |> x -> replace(x, re1 => "\$") |> x -> replace(x, re2 => s"$^\1$")
end

examples_list = ["symmetric_strip.jl", "resistive_square_patch.jl", 
                 "cross_on_dielectric_substrate.jl", 
                 "square_loop_absorber.jl", "band_pass_filter.jl",
                 "cpss_optimization.jl"]

# Adds examples as subsections to the Examples page:
function replace_includes(str)
  for ex in examples_list
      content = read(ex, String)
      str = replace(str, "include(\"$(ex)\")" => content)
  end
  return str
end

flist = ["manual.jl"]
for file in flist
    Literate.markdown(file, "../src", codefence=("```@repl manual" => "```"))
    Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
end

for file in examples_list
  Literate.markdown(file, ".", credit=false)
  Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
end

function postinclude(str)
  for file in examples_list
    mdfile = splitext(file)[1] * ".md"
    str *= read(mdfile, String)
  end
  str
end

Literate.markdown("examples.jl", "../src", postprocess=postinclude, credit=false)

