#!/usr/bin/env bash
# Requires sudo. System VTK9 on Ubuntu 24.04 uses Qt5 — install matching stack.
set -euo pipefail
sudo apt-get update
sudo apt-get install -y \
  qtbase5-dev qt5-qmake libqt5opengl5-dev \
  libvtk9-dev libvtk9-qt-dev \
  libopenmpi-dev \
  libgl1-mesa-dev libxkbcommon-dev \
  ninja-build cmake
echo "Done. Configure with Qt5 + system VTK:"
echo "  rm -rf software/desktop/build"
echo "  cmake -G Ninja -B software/desktop/build -S software/desktop -DCMAKE_BUILD_TYPE=Release"
echo "  cmake --build software/desktop/build -j\"\$(nproc)\""
