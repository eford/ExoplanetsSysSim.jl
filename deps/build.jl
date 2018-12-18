using Pkg, CORBITS

if !isfile(joinpath(dirname(pathof(CORBITS)),"../libcorbits.so"))
   Pkg.build("CORBITS")
end

if !Sys.iswindows()
   symlink( joinpath(dirname(pathof(CORBITS)),"../libcorbits.so"), joinpath(Pkg.devdir(),"ExoplanetsSysSim","libcorbits.so"))
else
   cp( joinpath(Pkg.devdir(),"CORBITS","libcorbits.so"), joinpath(Pkg.devdir(),"ExoplanetsSysSim","libcorbits.so"))
   cp( joinpath(dirname(pathof(CORBITS)),"../libcorbits.so"), joinpath(Pkg.devdir(),"ExoplanetsSysSim","libcorbits.so"))
end


