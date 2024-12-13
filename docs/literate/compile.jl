using Literate

function notebook_filter(str)
  re1 = r"(?<!`)``(?!`)"  # Two backquotes not preceded by nor followed by another
  re2 = r"\[\^(\d)\]:?" # e.g. "[^1]" or "[^1]:" while capturing the digit
  str |> x -> replace(x, re1 => "\$") |> x -> replace(x, re2 => s"$^\1$")
end

examples_list = ["symmetric_strip.jl", "resistive_square_patch.jl", 
                 "cross_on_dielectric_substrate.jl", 
                 "square_loop_absorber.jl", "flexible_absorber.jl", "splitringexample.jl", "reflectarray_example.jl",
                 "band_pass_filter.jl", "cpss1.jl", "cpss_optimization.jl", "cpss2.jl", "splitring_cpss.jl",
                 "angular_ss_example.jl", "tepfile_creation_example.jl"]


flist = ["manual.jl"]
for file in flist
    Literate.markdown(file, "../src", codefence=("```@repl manual" => "```"), credit=false)
    #Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
    Literate.notebook(file, "../notebooks", execute=false)
end

for (i,file) in enumerate(examples_list)
    fnpre = splitext(file)[1]
    Literate.markdown(file, "../src", credit=false)
    #Literate.notebook(file, "../notebooks", preprocess=notebook_filter, execute=false)
    Literate.notebook(file, "../notebooks", execute=false)
end

function postinclude(str)
  for file in examples_list
    mdfile = joinpath("..", "src", splitext(file)[1] * ".md")
    str *= read(mdfile, String)
    rm(mdfile)
  end
  return str
end

Literate.markdown("examples.jl", "../src", postprocess=postinclude, credit=false)
