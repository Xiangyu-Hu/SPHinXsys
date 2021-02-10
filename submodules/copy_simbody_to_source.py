import shutil

shutil.copytree('simbody/Simbody', '../SPHINXsys/src/shared/Simbody')
shutil.copytree('simbody/SimTKcommon', '../SPHINXsys/src/shared/SimTKcommon')
shutil.copytree('simbody/SimTKmath', '../SPHINXsys/src/shared/SimTKmath')

#find . -mindepth 0 -type d -exec cp ../CMakeLists.txt {} \;