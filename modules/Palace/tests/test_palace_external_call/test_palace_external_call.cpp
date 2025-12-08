#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
    std::cout << "[SPHinXsys-Palace] External Palace call test\n";

    // Palace 的主可执行文件路径，由 CMake 传进来
    const char *palace_exe = PALACE_EXECUTABLE;

    // 控制 Palace 的配置文件路径（在当前运行目录的 data/ 下面）
    std::string config_file = "data/rings.json";

    // 用 mpirun 调 Palace（你也可以改成 -np 4 之类的）
   std::string cmd = PALACE_EXECUTABLE;  // palace_exe
    cmd += " ";
    cmd += config_file;

    std::cout << "Running command:\n  " << cmd << std::endl;

    int ret = std::system(cmd.c_str());
    if (ret != 0)
    {
        std::cerr << "Palace process failed with code " << ret << std::endl;
        return ret;
    }

    std::cout << "Palace external run finished successfully.\n";
    return 0;
}
