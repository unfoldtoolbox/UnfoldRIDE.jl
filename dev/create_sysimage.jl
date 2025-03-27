# This script creates a sysimage for the RIDE project
# WARNING: This will take a LONG time (1-2 hours)
# The sysimage massively speeds up julia startup time
# After the sysimage is created, you need to add:
# "-JRIDE_Sysimage.so" to your "julia.additionalArgs" in your vscode settings
# You can use Base.loaded_modules to confirm that the sysimage is loaded correctly

include("./runner.jl")
# Get all packages from the current project
all_packages = String[]
for (uuid, dep) in Pkg.dependencies()
    dep.is_direct_dep || continue
    dep.version === nothing && continue
    push!(all_packages, dep.name)
end

# Remove unneeded packages
do_not_include = ["PackageCompiler"]
package_list = filter(x -> x ∉ do_not_include, all_packages)

using PackageCompiler

PackageCompiler.create_sysimage(
    package_list;
    sysimage_path = "RIDE_Sysimage.so",
    precompile_execution_file = "./dev/runner.jl",
)