project('hpdml', 'cpp',
	default_options : ['cpp_std=c++11','buildtype=release','prefix=/home/CAM/plewis/yubo'],
	version : '1.0')
cpp = meson.get_compiler('cpp')

# These lines specify the locations of the NCL, Boost, and BeagleLib compiled libraries
lib_system = cpp.find_library('boost_system', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_filesystem = cpp.find_library('boost_filesystem', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_program_options = cpp.find_library('boost_program_options', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_ncl = cpp.find_library('ncl', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_beagle = cpp.find_library('hmsbeagle', dirs: ['/home/CAM/plewis/lib'], required: true)

# These lines specify the locations of header files for the NCL, Boost, BeagleLib, and Eigen library
incl_beagle = include_directories('/home/CAM/plewis/include/libhmsbeagle-1')
incl_ncl = include_directories('/home/CAM/plewis/include')
incl_boost = include_directories('/home/CAM/plewis/boost_1_73_0')
incl_eigen = include_directories('/home/CAM/plewis/eigen-eigen-323c052e1731')

# This line creates the executable file
executable('hpdml', 'main.cpp', install: true, install_dir: '.', dependencies: [lib_beagle,lib_ncl,lib_program_options,lib_system,lib_filesystem], include_directories: [incl_beagle,incl_ncl,incl_boost,incl_eigen])

# These lines just copy files to the install directory
install_data('hpdml.conf', install_dir: '.')
install_data('rbcl10.nex', install_dir: '.')
install_data('rbcl10.tre', install_dir: '.')
install_data('s.sh', install_dir: '.')
install_data('test.sh', install_dir: '.')

