help([==[

Description
===========
MATLAB is a high-level language and interactive environment
 that enables you to perform computationally intensive tasks faster than with
 traditional programming languages such as C, C++, and Fortran.


More information
================
 - Homepage: https://www.mathworks.com/products/matlab
]==])

whatis([==[Description: MATLAB is a high-level language and interactive environment
 that enables you to perform computationally intensive tasks faster than with
 traditional programming languages such as C, C++, and Fortran.]==])
whatis([==[Homepage: https://www.mathworks.com/products/matlab]==])
whatis([==[URL: https://www.mathworks.com/products/matlab]==])

local root = "/opt/apps/testapps/el7/software/staging/MATLAB/2022a"

conflict("MATLAB")

prepend_path("CMAKE_PREFIX_PATH", root)
prepend_path("PATH", pathJoin(root, "bin"))
setenv("EBROOTMATLAB", root)
setenv("EBVERSIONMATLAB", "2022a")
setenv("EBDEVELMATLAB", pathJoin(root, "easybuild/MATLAB-2022a-easybuild-devel"))

setenv("LM_LICENSE_FILE", "48832@matlablm.shef.ac.uk")
prepend_path("PATH", root)
prepend_path("LD_LIBRARY_PATH", pathJoin(root, "runtime/glnxa64"))
prepend_path("LD_LIBRARY_PATH", pathJoin(root, "bin/glnxa64"))
prepend_path("LD_LIBRARY_PATH", pathJoin(root, "sys/os/glnxa64"))
setenv("_JAVA_OPTIONS", "-Xmx2048m")
-- Built with EasyBuild version 4.7.0
