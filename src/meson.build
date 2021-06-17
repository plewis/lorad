project('hpdml', 'cpp',
	default_options : ['cpp_std=c++11','buildtype=release','prefix=/home/CAM/plewis/yubo'],
	version : '1.0')
cpp = meson.get_compiler('cpp')
lib_regex = cpp.find_library('boost_regex', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_system = cpp.find_library('boost_system', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_filesystem = cpp.find_library('boost_filesystem', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_program_options = cpp.find_library('boost_program_options', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_ncl = cpp.find_library('ncl', dirs: ['/home/CAM/plewis/lib/static'], required: true)
lib_beagle = cpp.find_library('hmsbeagle', dirs: ['/home/CAM/plewis/lib'], required: true)
#lib_beagle = cpp.find_library('hmsbeagle', dirs: ['/home/CAM/plewis/lib/static'], required: true)
#lib_beagle_cpu = cpp.find_library('hmsbeagle-cpu', dirs: ['/home/CAM/plewis/lib/static'], required: true)
#lib_beagle_cpu_sse = cpp.find_library('hmsbeagle-cpu-sse', dirs: ['/home/CAM/plewis/lib/static'], required: true)
incl_beagle = include_directories('/home/CAM/plewis/include/libhmsbeagle-1')
incl_ncl = include_directories('/home/CAM/plewis/include')
incl_boost = include_directories('/home/CAM/plewis/boost_1_73_0')
incl_eigen = include_directories('/home/CAM/plewis/eigen-eigen-323c052e1731')
executable('hpdml', 'main.cpp', install: true, install_dir: '.', dependencies: [lib_beagle,lib_ncl,lib_program_options,lib_system,lib_filesystem,lib_regex], include_directories: [incl_beagle,incl_ncl,incl_boost,incl_eigen])
install_data('hpdml.conf', install_dir: '.')
install_data('rbcl10.nex', install_dir: '.')
install_data('rbcl10.tre', install_dir: '.')
#install_data('rbcl738.nex', install_dir: '.')
#install_data('rbcl738nj.tre', install_dir: '.')
#install_data('test.tre', install_dir: '.')
#install_data('rbcL.nex', install_dir: '.')
#install_data('rbcLjc.tre', install_dir: '.')
#install_data('go.sh', install_dir: '.')
#install_data('datatest.nex', install_dir: '.')


# c++ 
# -O3 
# -I$HOME/hpd-histogram/src 
# -I$HOME/boost_1_73_0 
# -I$HOME/include/libhmsbeagle-1 
# -I$HOME/include/
# -I$HOME/eigen-eigen-323c052e1731 
# -std=c++11
# -c main.cpp 
# $HOME/lib/static/libhmsbeagle.a 
# $HOME/lib/static/libhmsbeagle-cpu.a 
# $HOME/lib/static/libhmsbeagle-cpu-sse.a 
# $HOME/lib/static/libncl.a 
# $HOME/lib/static/libboost_program_options.a 
# $HOME/lib/static/libboost_system.a 
# $HOME/lib/static/libboost_filesystem.a 
# $HOME/lib/static/libboost_regex.a
# -o doof



