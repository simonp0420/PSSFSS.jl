"""
# Constants
Defines global constants used throughout the code, such as:
- speed of light
- permitivity of free space
- permeability of free space
- impedance of free space
- tolerances
- constants for ease of use
- various configuration constants 
"""
module Constants

export μ₀, ϵ₀, c₀, η₀, twopi, fourpi, tol, dbmin, tdigits

"Permeability of free space [H/m]"
const μ₀ = 1.25663706212 * 1e-6

"Speed of light [m/s]"
const c₀ = 299792458.0

"Permeability of free space [F/m]"
const ϵ₀ = 1 / (μ₀ * c₀^2)

"Intrinsic impedance of free space [F/m]"
const η₀ = sqrt(μ₀ / ϵ₀)

const twopi = 2π
const fourpi = 4π
const tol = 1e-4

"Minimum modal attenuation"
const dbmin = 30.0

"Min. elect. length for a layer to not be included in a `GBLOCK`"
const min_elength = 0.18 # Results in about 320 modes

"Number of digits to use in printing elapsed time"
const tdigits = 4

end
