for repo in opm-common opm-material opm-grid opm-models
do
    git clone https://github.com/daavid00/$repo.git --branch wettability
    mkdir $repo/build-cmake
    cd $repo/build-cmake
    cmake ..
    make -j 5
    cd ../..
done
cd opm-models/build-cmake
make -j 5 wa
cd ../..
git clone https://github.com/daavid00/opm-simulators.git --branch wettability
mkdir opm-simulators/build-cmake
cd opm-simulators/build-cmake
cmake ..
make -j 5 flow
cd ../..
git clone https://github.com/daavid00/kassa-et-al-2021/figures.git
