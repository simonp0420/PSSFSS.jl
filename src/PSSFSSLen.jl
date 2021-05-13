module PSSFSSLen
using Reexport

@reexport using Unitful
export mm, cm, inch, mil, PSSFSSLength


const mm = unit(1u"mm")
const cm = unit(1u"cm")
const inch = unit(1u"inch")
const mil = unit(1u"mil")
const PSSFSSLength = Union{typeof(mm),typeof(cm),typeof(inch),typeof(mil)}

end